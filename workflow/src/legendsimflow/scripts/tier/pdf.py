# ruff: noqa: I002

# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>,
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import argparse
import itertools
import logging
import re
from collections.abc import Mapping

import awkward as ak
import hist
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
import numpy as np
from lgdo import Histogram, Scalar, Struct
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, utils
from legendsimflow.metadata import get_sanitized_fccd_uncertainty, get_tier_settings
from legendsimflow.pdf import (
    aoe_categories,
    build_pdf_grid,
    category_levels,
    collected_energy,
    fccd_by_uid,
    fccd_categories,
    flatten_surface_steps,
    format_pdf_struct,
    grid_axis_values,
    lar_categories,
)
from legendsimflow.scripts import log_script_invocation


def _detector_group_uids(
    groups: Mapping[str, str] | None,
    name_to_uid: Mapping[str, int],
    log: logging.Logger,
) -> dict[str, np.ndarray]:
    """Return the UIDs of the detectors in each group, plus the implicit ``all`` group."""
    groups = dict(groups or {})
    if groups.get("all", ".*") != ".*":
        log.warning(
            "detector_groups sets 'all' to %r; it will be overridden by the implicit all-detector group '.*'",
            groups["all"],
        )
    groups["all"] = ".*"

    group_uids = {}
    for name, pattern in groups.items():
        rgx = re.compile(pattern)
        matched = sorted(uid for n, uid in name_to_uid.items() if rgx.fullmatch(n))
        if not matched:
            log.warning("detector group %r matches no detectors", name)
        group_uids[name] = np.asarray(matched, dtype=np.int64)
    return group_uids


def _fill_mul1(
    histograms: Mapping[tuple[str, ...], Mapping[str, hist.Hist]],
    energy: np.ndarray,
    energy_thr: float,
    hit_evt: np.ndarray,
    is_good: np.ndarray,
    in_category: Mapping[str, np.ndarray],
    in_group: Mapping[str, np.ndarray],
) -> None:
    """Fill the multiplicity-1 histograms with the HPGe hit energies `energy` (keV).

    All arrays have one entry per hit; `hit_evt` is the event of each hit.
    """
    # the multiplicity is counted here, not read from geds.multiplicity, so
    # that it matches the float32 or recomputed energies
    above_thr = energy > energy_thr
    multiplicity = np.bincount(hit_evt, weights=above_thr)
    n_not_good = np.bincount(hit_evt, weights=above_thr & ~is_good)

    # the one hit above threshold of the multiplicity-1 events, in an ON detector
    m1 = np.flatnonzero(
        above_thr & (multiplicity[hit_evt] == 1) & (n_not_good[hit_evt] == 0)
    )
    energy_m1 = energy[m1]
    for combo, per_group_hists in histograms.items():
        sel = np.ones(len(m1), dtype=bool)
        for c in combo:
            sel &= in_category[c][m1]
        for g, h in per_group_hists.items():
            h.fill(energy_m1[sel & in_group[g][m1]])


@snakemake_compatible(
    mapping={
        "cvt_file": "input.cvt_file",
        "pdf_file": "output[0]",
        "log_file": "log[0]",
        "simid": "wildcards.simid",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(description="Build the pdf tier.")
    parser.add_argument("--cvt-file", required=True, help="input cvt tier file")
    parser.add_argument("--pdf-file", required=True, help="output pdf tier file")
    parser.add_argument("--simid", required=True, help="simulation ID wildcard")
    parser.add_argument("--log-file", default=None, help="log file")
    parser.add_argument(
        "--simflow-config",
        "--config",
        dest="simflow_config",
        required=True,
        help="simflow config YAML path",
    )
    args = parser.parse_args()

    config = utils.init_simflow_context(args.simflow_config, workflow=None).config

    tier_pdf_settings = get_tier_settings(config, "pdf")
    tier_evt_settings = get_tier_settings(config, "evt")
    tier_hit_settings = get_tier_settings(config, "hit")

    cvt_file = str(nersc.dvs_ro(config, args.cvt_file))
    pdf_file = args.pdf_file

    log = ldfs.utils.build_log(config.metadata.simprod.config.logging, args.log_file)
    log_script_invocation(log, "tier-pdf", parser, args)

    # the number of simulated primary events is written by remage at the root
    # of each stp file, forwarded through the evt tier and summed over all jobs
    # at the cvt tier.
    number_of_primaries = int(lh5.read("number_of_simulated_events", cvt_file).value)
    log.info(
        "... read number of simulated events from %s: %s", cvt_file, number_of_primaries
    )

    output = Struct()

    # without HPGe data (skip_hit) there is nothing to histogram
    if not tier_evt_settings.skip_hit:
        has_spms_coinc = not tier_evt_settings.skip_opt
        # the A/E categories come from the pulse-shape library when it was
        # simulated, otherwise from the single template
        aoe_model = None
        if tier_hit_settings.simulate_psd_with_psl:
            aoe_model = "pulse_lib"
        elif tier_hit_settings.simulate_psd:
            aoe_model = "single_temp"
        fccd_thr = tier_pdf_settings.fccd_tag_energy_thr_kev
        if fccd_thr is not None and not tier_hit_settings.store_surface_steps:
            msg = (
                "fccd_tag_energy_thr_kev needs the steps near the n+ surface, "
                "run the hit tier with store_surface_steps"
            )
            raise ValueError(msg)
        morph = tier_pdf_settings.morph_dead_layer
        if morph and not tier_hit_settings.store_surface_steps:
            msg = (
                "morph_dead_layer needs the steps near the n+ surface, "
                "run the hit tier with store_surface_steps"
            )
            raise ValueError(msg)

        # name to UID map written at evt-tier build time
        detector_uids = lh5.read("detector_uids", cvt_file)
        name_to_uid = {n: int(detector_uids[n].value) for n in detector_uids}
        group_uids = _detector_group_uids(
            tier_pdf_settings.detector_groups, name_to_uid, log
        )
        if fccd_thr is not None or morph:
            sorted_uids, fccd_mm = fccd_by_uid(
                config.metadata,
                name_to_uid,
                tier_hit_settings.surface_depth_in_mm,
            )

        # one histogram per combination of categories, one per level, and
        # detector group; the combinations do not overlap
        levels = category_levels(has_spms_coinc, aoe_model, fccd_thr)
        category_desc = {k: v for level in levels for k, v in level.items()}
        combos = list(itertools.product(*levels))
        paths = {combo: "/".join(("mul1", *combo)) for combo in combos}
        descriptions = {
            paths[combo]: "multiplicity 1"
            + "".join(f", {category_desc[c]}" for c in combo)
            for combo in combos
        }

        # 1 keV/bin over 0-6000 keV; boost-histogram only provides Double
        # (float64) storage
        histograms = {
            combo: {g: hist.new.Reg(6000, 0, 6000).Double() for g in group_uids}
            for combo in combos
        }
        # 2-D: (E_min, E_max) for multiplicity-2 events; not split by group
        mul2_hist = hist.new.Reg(6000, 0, 6000).Reg(6000, 0, 6000).Double()

        # multiplicity-1 histograms at each point of the grid of dead-layer
        # models: the FCCD is shifted by fccd_sigma times its (negative or
        # positive) uncertainty, limited to the physical range
        if morph:
            uid_to_name = {uid: n for n, uid in name_to_uid.items()}
            fccd_unc_mm = np.array(
                [
                    get_sanitized_fccd_uncertainty(config.metadata, uid_to_name[u])
                    for u in sorted_uids
                ]
            )
            morph_axes = {
                name: tuple(tier_pdf_settings.morph_grid[name])
                for name in ("fccd_sigma", "dlf")
            }
            morph_fccd_range = (
                tier_pdf_settings.morph_min_fccd_in_mm,
                tier_hit_settings.surface_depth_in_mm,
            )
            morph_points = {}
            for sigma, dlf in itertools.product(
                *(grid_axis_values(*axis) for axis in morph_axes.values())
            ):
                unc = fccd_unc_mm[:, 1] if sigma >= 0 else fccd_unc_mm[:, 0]
                fccd_morph = np.clip(fccd_mm + sigma * unc, *morph_fccd_range)
                morph_points[(float(sigma), float(dlf))] = (
                    fccd_morph,
                    {
                        combo: {
                            g: hist.new.Reg(6000, 0, 6000).Double() for g in group_uids
                        }
                        for combo in combos
                    },
                )
            log.info("... building the PDFs at %d dead-layer models", len(morph_points))

        log.info("... beginning iteration over cvt file")

        # read only the fields used here: the SiPM table and the surface steps
        # are large
        fields = ["geds/energy", "geds/rawid", "geds/quality/is_good_channel"]
        if has_spms_coinc:
            fields.append("coincident/spms")
        if aoe_model is not None:
            fields += ["geds/psd/is_good", f"geds/psd/{aoe_model}"]
        if fccd_thr is not None or morph:
            fields += [
                "geds/surface/edep",
                "geds/surface/distance_to_nplus",
                "geds/surface/energy_bulk",
            ]

        # a single cvt file might be large, need to iterate
        for chunk in lh5.LH5Iterator(
            cvt_file,
            "evt",
            field_mask=fields,
            buffer_len=tier_pdf_settings.buffer_len,
        ):
            data = chunk.view_as("ak")
            geds = data.geds

            # apply the energy threshold again
            above_thr = geds.energy > tier_evt_settings.geds_energy_thr_kev
            multiplicity = ak.sum(above_thr, axis=-1)
            is_good = ak.all(geds.quality.is_good_channel[above_thr], axis=-1)

            energy_m2 = geds.energy[above_thr][(multiplicity == 2) & is_good]
            mul2_hist.fill(ak.min(energy_m2, axis=-1), ak.max(energy_m2, axis=-1))

            # one entry per hit, the categories and groups do not depend on
            # the dead-layer model
            hit_evt = np.repeat(np.arange(len(data)), ak.to_numpy(ak.num(geds.energy)))
            energy = ak.to_numpy(ak.flatten(geds.energy))
            uid = ak.to_numpy(ak.flatten(geds.rawid))
            is_good = ak.to_numpy(ak.flatten(geds.quality.is_good_channel))
            in_category = {}
            if has_spms_coinc:
                in_category |= {k: v[hit_evt] for k, v in lar_categories(data).items()}
            if aoe_model is not None:
                in_category |= aoe_categories(geds.psd, aoe_model)
            in_group = {g: np.isin(uid, uids) for g, uids in group_uids.items()}
            fill_args = (
                tier_evt_settings.geds_energy_thr_kev,
                hit_evt,
                is_good,
            )

            if fccd_thr is not None or morph:
                steps = flatten_surface_steps(geds, sorted_uids)
            if fccd_thr is not None:
                fccd_tag = fccd_categories(steps, fccd_mm, fccd_thr)
            _fill_mul1(
                histograms,
                energy,
                *fill_args,
                in_category | fccd_tag if fccd_thr is not None else in_category,
                in_group,
            )

            if morph:
                collected_nominal = collected_energy(
                    steps, fccd_mm, tier_hit_settings.dead_layer_fraction
                )
                for (_, dlf), (fccd_morph, hists) in morph_points.items():
                    # nominal (smeared) energy plus the change of collected energy
                    energy_morph = energy + (
                        collected_energy(steps, fccd_morph, dlf) - collected_nominal
                    )
                    if fccd_thr is not None:
                        fccd_tag = fccd_categories(steps, fccd_morph, fccd_thr)
                    _fill_mul1(
                        hists,
                        energy_morph,
                        *fill_args,
                        in_category | fccd_tag if fccd_thr is not None else in_category,
                        in_group,
                    )

        log.info("... convert histograms to lgdo")
        output = format_pdf_struct(
            {paths[combo]: h for combo, h in histograms.items()}, descriptions
        )
        output.add_field(
            "mul2/all",
            Histogram(
                mul2_hist,
                attrs={
                    "description": "multiplicity-2 events: (E_low, E_high) for pairs of ON detectors"
                },
            ),
        )

    lh5.write(output, "pdf", pdf_file, wo_mode="w")
    if not tier_evt_settings.skip_hit and morph:
        log.info("... convert dead-layer grids to lgdo")
        grids = build_pdf_grid(
            lambda pars: {
                f"{paths[combo]}/{g}": Histogram(
                    h, attrs={"description": descriptions[paths[combo]]}
                )
                for combo, per_group in morph_points[(pars["fccd_sigma"], pars["dlf"])][
                    1
                ].items()
                for g, h in per_group.items()
            },
            morph_axes,
        )
        morph_output = Struct()
        for path, grid in grids.items():
            morph_output.add_field(path, grid)
        lh5.write(morph_output, "pdf_morph", pdf_file, wo_mode="append")
    lh5.write(Scalar(number_of_primaries), "nr_sim_events", pdf_file, wo_mode="append")


if __name__ == "__main__":
    main()
