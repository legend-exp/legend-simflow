# ruff: noqa: I002

# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>
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

import awkward as ak
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
from lgdo import Scalar
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, utils
from legendsimflow.metadata import get_tier_settings, runinfo
from legendsimflow.pdf import (
    fill_histograms,
    format_pdf_struct,
    group_detectors,
    new_histograms,
)
from legendsimflow.scripts import log_script_invocation
from legendsimflow.selections import (
    aoe_categories,
    aoe_model,
    category_levels,
    data_aoe_flags,
    data_lar_veto,
    is_clean_data_event,
    lar_categories,
    one_hit,
    select_hits,
)


@snakemake_compatible(
    mapping={
        "runid": "wildcards.runid",
        "evt_files": "input.evt_files",
        "modeled_detectors": "params.modeled_detectors",
        "pdf_file": "output[0]",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Histogram the LEGEND-200 data of a run in the format of the pdf tier."
    )
    parser.add_argument(
        "--runid", required=True, help="run ID (e.g. l200-p03-r000-phy)"
    )
    parser.add_argument(
        "--evt-files",
        nargs="+",
        default=None,
        help="evt (or pet) tier files of the run, default: those in paths.l200data",
    )
    parser.add_argument(
        "--modeled-detectors",
        nargs="*",
        default=[],
        help="names of the HPGe detectors whose A/E the simulation models",
    )
    parser.add_argument("--pdf-file", required=True, help="output file")
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

    log = ldfs.utils.build_log(config.metadata.simprod.config.logging, args.log_file)
    log_script_invocation(log, "tier-pdf-l200data", parser, args)

    evt_files = args.evt_files
    if evt_files is None:
        l200data = config.paths.l200data
        evt_files = utils.lookup_evt_files(
            l200data, args.runid, utils.get_evt_tier_name(l200data)
        )
    evt_files = [str(nersc.dvs_ro(config, f)) for f in evt_files]

    meta = config.metadata
    info = runinfo(meta, args.runid)
    name_to_rawid = {
        name: int(ch.daq.rawid)
        for name, ch in meta.channelmap(info.start_key).items()
        if ch.get("system") == "geds"
    }
    group_rawids = group_detectors(tier_pdf_settings.detector_groups, name_to_rawid)
    modeled_rawids = [name_to_rawid[n] for n in args.modeled_detectors]

    has_spms_coinc = not tier_evt_settings.skip_opt
    model = aoe_model(tier_hit_settings)
    levels = category_levels(has_spms_coinc, model, None)
    histograms, mul2_hist = new_histograms(levels, group_rawids)

    fields = [
        "trigger/is_forced",
        "coincident/puls",
        "coincident/muon",
        "coincident/muon_offline",
        "geds/energy",
        "geds/rawid",
        "geds/quality/is_good_channel",
        "geds/quality/is_bb_like",
    ]
    if has_spms_coinc:
        fields += ["spms/multiplicity", "spms/energy_sum"]
    if model is not None:
        fields += ["geds/psd/is_good", "geds/psd/low_aoe/is_single_site"]
    if model == "pulse_lib":
        fields.append("geds/psd/high_aoe/is_bulk")

    log.info("... beginning iteration over %s", evt_files)

    for chunk in lh5.LH5Iterator(
        evt_files,
        "evt",
        field_mask=fields,
        buffer_len=tier_pdf_settings.buffer_len,
    ):
        data = chunk.view_as("ak")
        data = data[is_clean_data_event(data)]
        is_hit, multiplicity, is_good = select_hits(
            data, tier_evt_settings.geds_energy_thr_kev
        )

        energy_m2 = data.geds.energy[is_hit][(multiplicity == 2) & is_good]
        mul2_hist.fill(ak.min(energy_m2, axis=-1), ak.max(energy_m2, axis=-1))

        is_m1 = (multiplicity == 1) & is_good
        evt_m1 = data[is_m1]
        hit_m1 = is_hit[is_m1]
        energy = ak.to_numpy(one_hit(evt_m1.geds.energy, hit_m1))
        rawid = ak.to_numpy(one_hit(evt_m1.geds.rawid, hit_m1))

        in_category = {}
        if has_spms_coinc:
            # the flag is not read from the data production: it must follow
            # the definition of the Simflow evt tier
            in_category |= lar_categories(
                data_lar_veto(
                    evt_m1,
                    tier_evt_settings.lar_veto_multiplicity_thr,
                    tier_evt_settings.lar_veto_energy_sum_pe_thr,
                )
            )
        if model is not None:
            in_category |= aoe_categories(
                *data_aoe_flags(evt_m1, hit_m1, model, modeled_rawids)
            )

        fill_histograms(histograms, energy, rawid, in_category, group_rawids)

    log.info("... convert histograms to lgdo")
    output = format_pdf_struct(histograms, levels, mul2_hist)

    lh5.write(output, "pdf", args.pdf_file, wo_mode="w")
    lh5.write(
        Scalar(float(info.livetime_in_s)),
        "livetime_in_s",
        args.pdf_file,
        wo_mode="append",
    )


if __name__ == "__main__":
    main()
