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
import copy
import itertools

import awkward as ak
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
import numpy as np
from lgdo import Scalar, Struct
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, utils
from legendsimflow.metadata import get_sanitized_fccd_uncertainty, get_tier_settings
from legendsimflow.pdf import (
    build_pdf_grid,
    collected_energy,
    fccd_layer_energy,
    fill_histograms,
    flatten_surface_steps,
    format_pdf_struct,
    grid_axis_values,
    group_detectors,
    mul1_hits,
    new_histograms,
)
from legendsimflow.scripts import log_script_invocation
from legendsimflow.selections import (
    aoe_categories,
    aoe_flags,
    aoe_model,
    category_levels,
    fccd_by_uid,
    fccd_categories,
    has_hpge_multiplicity,
    is_hpge_hit_above_thr,
    is_spms_coincident,
    lar_categories,
    one_hit,
)


def _aoe_flags_per_hit(
    psd: ak.Array, model: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Return the A/E flags of :func:`.selections.aoe_flags` for every simulated HPGe hit."""

    def flat(field):
        return ak.to_numpy(ak.flatten(field))

    aoe = psd[model]
    high = flat(aoe.is_high_aoe) if model == "pulse_lib" else None
    return flat(psd.is_good & aoe.has_aoe), flat(aoe.is_single_site), high


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

    log = ldfs.utils.build_log(
        config.metadata.simprod.config.logging.to_dict(), args.log_file
    )
    log_script_invocation(log, "tier-pdf", parser, args)

    # the number of simulated primary events is written by remage at the root
    # of each stp file, forwarded through the evt tier and summed over all jobs
    # at the cvt tier.
    number_of_primaries = int(lh5.read("number_of_simulated_events", cvt_file).value)
    log.info(
        "... read number of simulated events from %s: %s", cvt_file, number_of_primaries
    )

    output = Struct()
    morph_output = None

    # without HPGe data (skip_hit) there is nothing to histogram
    if not tier_evt_settings.skip_hit:
        has_spms_coinc = not tier_evt_settings.skip_opt
        model = aoe_model(tier_hit_settings)
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
        group_uids = group_detectors(tier_pdf_settings.detector_groups, name_to_uid)
        if fccd_thr is not None or morph:
            sorted_uids, fccd_mm = fccd_by_uid(
                config.metadata,
                name_to_uid,
                tier_hit_settings.surface_depth_in_mm,
            )

        levels = category_levels(has_spms_coinc, model, fccd_thr)
        histograms, mul2_hist = new_histograms(levels, group_uids)

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
                morph_points[(float(sigma), float(dlf))] = (
                    np.clip(fccd_mm + sigma * unc, *morph_fccd_range),
                    copy.deepcopy(histograms),
                )
            log.info("... building the PDFs at %d dead-layer models", len(morph_points))

        log.info("... beginning iteration over cvt file")

        # read only the fields used here: the SiPM table and the surface steps
        # are large
        fields = ["geds/energy", "geds/rawid", "geds/quality/is_good_channel"]
        if has_spms_coinc:
            fields += ["spms/multiplicity", "spms/energy_sum"]
        if model is not None:
            fields += ["geds/psd/is_good", f"geds/psd/{model}"]
        if fccd_thr is not None or morph:
            fields += ["geds/surface/edep", "geds/surface/distance_to_nplus"]
        if morph:
            fields.append("geds/surface/energy_bulk")

        # a single cvt file might be large, need to iterate
        for chunk in lh5.LH5Iterator(
            cvt_file,
            "evt",
            field_mask=fields,
            buffer_len=tier_pdf_settings.buffer_len,
        ):
            data = chunk.view_as("ak")
            thr = tier_evt_settings.geds_energy_thr_kev
            is_hit = is_hpge_hit_above_thr(data, thr)

            energy_m2 = data.geds.energy[is_hit][
                has_hpge_multiplicity(data, 2, None, thr)
            ]
            mul2_hist.fill(ak.min(energy_m2, axis=-1), ak.max(energy_m2, axis=-1))

            # multiplicity-1 events have exactly one hit above threshold
            evt_m1 = data[has_hpge_multiplicity(data, 1, None, thr)]
            hit = one_hit(evt_m1, thr)
            energy = ak.to_numpy(hit.energy)
            uid = ak.to_numpy(hit.rawid)

            # event masks of the categories of all levels
            in_category = {}
            if has_spms_coinc:
                in_category |= lar_categories(
                    is_spms_coincident(
                        evt_m1,
                        tier_evt_settings.lar_veto_multiplicity_thr,
                        tier_evt_settings.lar_veto_energy_sum_pe_thr,
                    )
                )
            if model is not None:
                in_category |= aoe_categories(*aoe_flags(evt_m1, model, thr))
            if fccd_thr is not None:
                fccd_of_hit = fccd_mm[np.searchsorted(sorted_uids, uid)]
                in_category |= fccd_categories(evt_m1, fccd_of_hit, fccd_thr, thr)

            fill_histograms(histograms, energy, uid, in_category, group_uids)

            if morph:
                # one entry per hit: the threshold, the multiplicity and the
                # FCCD tag change with the dead-layer model, the LAr and A/E
                # categories do not
                geds = data.geds
                hit_evt = np.repeat(
                    np.arange(len(data)), ak.to_numpy(ak.num(geds.energy))
                )
                energy_hit = ak.to_numpy(ak.flatten(geds.energy))
                uid_hit = ak.to_numpy(ak.flatten(geds.rawid))
                is_good = ak.to_numpy(ak.flatten(geds.quality.is_good_channel))
                in_category_hit = {}
                if has_spms_coinc:
                    lar = lar_categories(
                        is_spms_coincident(
                            data,
                            tier_evt_settings.lar_veto_multiplicity_thr,
                            tier_evt_settings.lar_veto_energy_sum_pe_thr,
                        )
                    )
                    in_category_hit |= {k: v[hit_evt] for k, v in lar.items()}
                if model is not None:
                    in_category_hit |= aoe_categories(
                        *_aoe_flags_per_hit(geds.psd, model)
                    )

                steps = flatten_surface_steps(geds, sorted_uids)
                collected_nominal = collected_energy(
                    steps, fccd_mm, tier_hit_settings.dead_layer_fraction
                )
                for (_, dlf), (fccd_morph, hists) in morph_points.items():
                    # nominal (smeared) energy plus the change of collected energy
                    energy_morph = energy_hit + (
                        collected_energy(steps, fccd_morph, dlf) - collected_nominal
                    )
                    in_cat = in_category_hit
                    if fccd_thr is not None:
                        tagged = fccd_layer_energy(steps, fccd_morph) > fccd_thr
                        in_cat = in_cat | {
                            "is_fccd_tagged": tagged,
                            "not_fccd_tagged": ~tagged,
                        }
                    m1 = mul1_hits(energy_morph, hit_evt, is_good, thr)
                    fill_histograms(
                        hists,
                        energy_morph[m1],
                        uid_hit[m1],
                        {k: v[m1] for k, v in in_cat.items()},
                        group_uids,
                    )

        log.info("... convert histograms to lgdo")
        output = format_pdf_struct(histograms, levels, mul2_hist)

        if morph:
            log.info("... convert dead-layer grids to lgdo")
            paths = [
                "/".join(("mul1", *combo, g))
                for combo, per_group in histograms.items()
                for g in per_group
            ]

            def morph_pdfs(pars):
                hists = morph_points[(pars["fccd_sigma"], pars["dlf"])][1]
                point = format_pdf_struct(hists, levels)
                return {p: point[p] for p in paths}

            morph_output = Struct()
            for path, grid in build_pdf_grid(morph_pdfs, morph_axes).items():
                morph_output.add_field(path, grid)

    lh5.write(output, "pdf", pdf_file, wo_mode="w")
    if morph_output is not None:
        lh5.write(morph_output, "pdf_morph", pdf_file, wo_mode="append")
    lh5.write(Scalar(number_of_primaries), "nr_sim_events", pdf_file, wo_mode="append")


if __name__ == "__main__":
    main()
