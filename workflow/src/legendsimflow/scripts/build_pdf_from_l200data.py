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
from pathlib import Path

import awkward as ak
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
from dbetto.utils import load_dict
from lgdo import Scalar
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, patterns, utils
from legendsimflow.exceptions import SimflowConfigError
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
    aoe_flags,
    aoe_model,
    attach_modelable_flag,
    category_levels,
    has_hpge_multiplicity,
    is_clean_event,
    is_hpge_hit_above_thr,
    is_spms_coincident,
    lar_categories,
    one_hit,
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
        default=None,
        help="names of the HPGe detectors whose A/E the simulation models, "
        "default: those of the run in the modelable-HPGe cache (pars/detinfo)",
    )
    parser.add_argument(
        "--pdf-file",
        default=None,
        help="output file, default: tier/pdf/l200data/<runid>-tier_pdf.lh5",
    )
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

    pdf_file = args.pdf_file
    if pdf_file is None:
        pdf_file = patterns.output_l200data_pdf_filename(config, runid=args.runid)
        Path(pdf_file).parent.mkdir(parents=True, exist_ok=True)

    # the data production stores only the hits above 25 keV
    energy_thr_kev = tier_evt_settings.geds_energy_thr_kev
    if energy_thr_kev < 25:
        msg = f"geds_energy_thr_kev is {energy_thr_kev}, the data has no hits below 25 keV"
        raise SimflowConfigError(msg)

    if args.modeled_detectors is None:
        cache = patterns.detinfo_filename(config, "is_modelable")
        is_modelable = load_dict(nersc.dvs_ro(config, cache))
        if args.runid not in is_modelable:
            msg = f"run {args.runid} is not in the modelable-HPGe cache {cache}"
            raise SimflowConfigError(msg)
    else:
        is_modelable = {args.runid: dict.fromkeys(args.modeled_detectors, True)}
    log.info(
        "... %d detectors with a modeled A/E", sum(is_modelable[args.runid].values())
    )

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

    has_spms_coinc = not tier_evt_settings.skip_opt
    model = aoe_model(tier_hit_settings)
    levels = category_levels(has_spms_coinc, model, None)
    histograms, mul2_hist = new_histograms(levels, group_rawids)

    fields = [
        "trigger/period",
        "trigger/run",
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
        data = data[is_clean_event(data)]
        if model is not None:
            data = attach_modelable_flag(data, config, is_modelable)
        is_hit = is_hpge_hit_above_thr(data, energy_thr_kev)

        energy_m2 = data.geds.energy[is_hit][
            has_hpge_multiplicity(data, 2, None, energy_thr_kev)
        ]
        mul2_hist.fill(ak.min(energy_m2, axis=-1), ak.max(energy_m2, axis=-1))

        evt_m1 = data[has_hpge_multiplicity(data, 1, None, energy_thr_kev)]
        hit = one_hit(evt_m1, energy_thr_kev)
        energy = ak.to_numpy(hit.energy)
        rawid = ak.to_numpy(hit.rawid)

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
            in_category |= aoe_categories(*aoe_flags(evt_m1, model, energy_thr_kev))

        fill_histograms(histograms, energy, rawid, in_category, group_rawids)

    log.info("... convert histograms to lgdo")
    output = format_pdf_struct(histograms, levels, mul2_hist)

    lh5.write(output, "pdf", pdf_file, wo_mode="w")
    lh5.write(
        Scalar(float(info.livetime_in_s)), "livetime_in_s", pdf_file, wo_mode="append"
    )


if __name__ == "__main__":
    main()
