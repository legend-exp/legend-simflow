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
from legendsimflow.metadata import get_tier_settings
from legendsimflow.pdf import (
    aoe_categories,
    category_levels,
    fccd_by_uid,
    fccd_categories,
    format_pdf_struct,
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

        # name to UID map written at evt-tier build time
        detector_uids = lh5.read("detector_uids", cvt_file)
        name_to_uid = {n: int(detector_uids[n].value) for n in detector_uids}
        group_uids = _detector_group_uids(
            tier_pdf_settings.detector_groups, name_to_uid, log
        )
        if fccd_thr is not None:
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

        # 1 keV/bin over 0-6000 keV; boost-histogram only provides Double
        # (float64) storage
        histograms = {
            combo: {g: hist.new.Reg(6000, 0, 6000).Double() for g in group_uids}
            for combo in combos
        }
        # 2-D: (E_min, E_max) for multiplicity-2 events; not split by group
        mul2_hist = hist.new.Reg(6000, 0, 6000).Reg(6000, 0, 6000).Double()

        log.info("... beginning iteration over cvt file")

        # read only the fields used here: the SiPM table and the surface steps
        # are large
        fields = ["geds/energy", "geds/rawid", "geds/quality/is_good_channel"]
        if has_spms_coinc:
            fields.append("coincident/spms")
        if aoe_model is not None:
            fields += ["geds/psd/is_good", f"geds/psd/{aoe_model}"]
        if fccd_thr is not None:
            fields += ["geds/surface/edep", "geds/surface/distance_to_nplus"]

        # a single cvt file might be large, need to iterate
        for chunk in lh5.LH5Iterator(
            cvt_file,
            "evt",
            field_mask=fields,
            buffer_len=tier_pdf_settings.buffer_len,
        ):
            data = chunk.view_as("ak")
            geds = data.geds

            # apply the energy threshold again. The multiplicity is counted
            # here, not read from geds.multiplicity, so that it matches the
            # float32 energies.
            above_thr = geds.energy > tier_evt_settings.geds_energy_thr_kev
            multiplicity = ak.sum(above_thr, axis=-1)
            is_good = ak.all(geds.quality.is_good_channel[above_thr], axis=-1)

            energy_m2 = geds.energy[above_thr][(multiplicity == 2) & is_good]
            mul2_hist.fill(ak.min(energy_m2, axis=-1), ak.max(energy_m2, axis=-1))

            # multiplicity-1 events have exactly one hit above threshold
            is_m1 = (multiplicity == 1) & is_good
            evt_m1 = data[is_m1]
            hit_m1 = above_thr[is_m1]
            energy = ak.to_numpy(evt_m1.geds.energy[hit_m1][:, 0])
            uid = ak.to_numpy(evt_m1.geds.rawid[hit_m1][:, 0])

            # event masks of the categories of all levels
            in_category = {}
            if has_spms_coinc:
                in_category |= lar_categories(evt_m1)
            if aoe_model is not None:
                in_category |= aoe_categories(evt_m1, hit_m1, aoe_model)
            if fccd_thr is not None:
                fccd_of_hit = fccd_mm[np.searchsorted(sorted_uids, uid)]
                in_category |= fccd_categories(evt_m1, hit_m1, fccd_of_hit, fccd_thr)

            in_group = {g: np.isin(uid, uids) for g, uids in group_uids.items()}
            for combo, per_group_hists in histograms.items():
                sel = np.ones(len(energy), dtype=bool)
                for c in combo:
                    sel &= in_category[c]
                for g, h in per_group_hists.items():
                    h.fill(energy[sel & in_group[g]])

        log.info("... convert histograms to lgdo")
        output = format_pdf_struct(
            {paths[combo]: h for combo, h in histograms.items()},
            {
                paths[combo]: "multiplicity 1"
                + "".join(f", {category_desc[c]}" for c in combo)
                for combo in combos
            },
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
    lh5.write(Scalar(number_of_primaries), "nr_sim_events", pdf_file, wo_mode="append")


if __name__ == "__main__":
    main()
