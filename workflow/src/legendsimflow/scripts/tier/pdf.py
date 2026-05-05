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
import re

import awkward as ak
import hist
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import numpy as np
from lgdo import Histogram, Scalar, Struct, lh5
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, utils
from legendsimflow.metadata import get_simconfig, get_tier_settings
from legendsimflow.scripts import log_script_invocation


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
    buffer_len = tier_pdf_settings.buffer_len

    cvt_file = nersc.dvs_ro(config, args.cvt_file)
    pdf_file = args.pdf_file
    log_file = args.log_file
    simid = args.simid

    log = ldfs.utils.build_log(config.metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "tier-pdf", parser, args)

    log.info("... reading data from %s", cvt_file)

    # detect what data is available in the cvt file
    _cvt_file_str = str(cvt_file)
    _evt_groups = {k.removeprefix("evt/") for k in lh5.ls(_cvt_file_str, "evt/")}
    has_geds = "geds" in _evt_groups
    has_spms_coinc = "coincident" in _evt_groups and "spms" in {
        k.removeprefix("evt/coincident/")
        for k in lh5.ls(_cvt_file_str, "evt/coincident/")
    }

    # resolve detector-group regexes against the name->uid map written at
    # evt-tier build time; the implicit "all" group is always emitted.
    detector_uids = lh5.read("detector_uids", str(cvt_file))
    name_to_uid = {name: int(detector_uids[name].value) for name in detector_uids}

    groups: dict[str, str] = dict(tier_pdf_settings.get("detector_groups", {}))
    if "all" in groups:
        log.warning(
            "detector_groups contains an 'all' key; it will be overridden by the implicit all-detector group"
        )
    groups["all"] = ".*"

    group_uids: dict[str, np.ndarray] = {}
    for gname, pat in groups.items():
        rgx = re.compile(pat)
        matched = sorted(uid for n, uid in name_to_uid.items() if rgx.fullmatch(n))
        # skip_hit upstream produces an empty detector_uids by design; warning would be spurious.
        if not matched and has_geds:
            log.warning("detector group %r matches no detectors", gname)
        group_uids[gname] = np.asarray(matched, dtype=np.int64)

    # a single cvt file might be large, need to iterate
    iterator = lh5.LH5Iterator(
        str(cvt_file),
        "evt",
        buffer_len=buffer_len,
    )

    # TODO: in future, read number_of_primaries directly from the stp file
    simconfig_block = get_simconfig(config, "stp", simid)
    number_of_primaries = (
        simconfig_block.primaries_per_job * simconfig_block.number_of_jobs
    )

    log.info("... extracted number of primaries from simconfig %s", number_of_primaries)

    # 1 keV/bin over 0-6000 keV; boost-histogram only provides Double (float64) storage
    h1 = hist.new.Reg(6000, 0, 6000).Double

    # 1-D histograms are split per detector group; mul2 stays global.
    histograms: dict[str, dict[str, hist.Hist]] = {}
    fail_histograms: dict[str, dict[str, hist.Hist]] = {}
    mul2_hist = None

    if has_geds:
        for cut in ("hit", "mul", "mul_psd"):
            histograms[cut] = {g: h1() for g in groups}
        fail_histograms["psd"] = {g: h1() for g in groups}
        # 2-D: (E_min, E_max) for multiplicity-2 events; not split by group
        mul2_hist = hist.new.Reg(6000, 0, 6000).Reg(6000, 0, 6000).Double()

    if has_spms_coinc:
        for cut in ("mul_lar", "mul_lar_psd"):
            histograms[cut] = {g: h1() for g in groups}
        fail_histograms["lar"] = {g: h1() for g in groups}

    log.info("... beginning iteration over cvt file")

    def _fill_per_group(
        per_group_hists: dict[str, hist.Hist],
        energy: ak.Array,
        det_mask: dict[str, ak.Array],
    ) -> None:
        for g in groups:
            per_group_hists[g].fill(ak.flatten(energy[det_mask[g]]))

    # events without a valid or simulated A/E are classified as background (cut).
    def _psd_mask(d: ak.Array) -> ak.Array:
        return ak.all(
            d.geds.psd.is_good & d.geds.psd.has_aoe & d.geds.psd.is_single_site,
            axis=-1,
        )

    for chunk in iterator:
        data = chunk.view_as("ak")

        if has_geds:
            # awkward 2.x has no ak.isin: build per-group masks via flatten / np.isin / unflatten.
            counts = ak.num(data.geds.rawid)
            flat_rawid = ak.flatten(data.geds.rawid)
            det_in_group = {
                g: ak.unflatten(np.isin(flat_rawid, uids), counts)
                for g, uids in group_uids.items()
            }

            good_channel_mask = ak.all(data.geds.is_good_channel, axis=-1)

            data_hit = data[good_channel_mask]
            det_hit = {g: m[good_channel_mask] for g, m in det_in_group.items()}
            _fill_per_group(histograms["hit"], data_hit.geds.energy, det_hit)

            m1_event_mask = (data.geds.multiplicity == 1) & good_channel_mask
            data_m1 = data[m1_event_mask]
            det_m1 = {g: m[m1_event_mask] for g, m in det_in_group.items()}
            _fill_per_group(histograms["mul"], data_m1.geds.energy, det_m1)

            data_m1_lar = data_m1
            det_m1_lar = det_m1
            if has_spms_coinc:
                lar_pass_m1 = ~data_m1.coincident.spms
                data_m1_lar = data_m1[lar_pass_m1]
                det_m1_lar = {g: m[lar_pass_m1] for g, m in det_m1.items()}
                _fill_per_group(
                    histograms["mul_lar"], data_m1_lar.geds.energy, det_m1_lar
                )

                lar_fail_m1 = data_m1.coincident.spms
                data_m1_lar_fail = data_m1[lar_fail_m1]
                det_m1_lar_fail = {g: m[lar_fail_m1] for g, m in det_m1.items()}
                _fill_per_group(
                    fail_histograms["lar"],
                    data_m1_lar_fail.geds.energy,
                    det_m1_lar_fail,
                )

            psd_pass_m1 = _psd_mask(data_m1)
            data_m1_psd = data_m1[psd_pass_m1]
            det_m1_psd = {g: m[psd_pass_m1] for g, m in det_m1.items()}
            _fill_per_group(histograms["mul_psd"], data_m1_psd.geds.energy, det_m1_psd)

            if has_spms_coinc:
                psd_pass_m1_lar = _psd_mask(data_m1_lar)
                data_m1_lar_psd = data_m1_lar[psd_pass_m1_lar]
                det_m1_lar_psd = {g: m[psd_pass_m1_lar] for g, m in det_m1_lar.items()}
                _fill_per_group(
                    histograms["mul_lar_psd"],
                    data_m1_lar_psd.geds.energy,
                    det_m1_lar_psd,
                )

            psd_fail_mask = (
                ak.all(data_m1.geds.psd.is_good & data_m1.geds.psd.has_aoe, axis=-1)
                & ~psd_pass_m1
            )
            data_m1_psd_fail = data_m1[psd_fail_mask]
            det_m1_psd_fail = {g: m[psd_fail_mask] for g, m in det_m1.items()}
            _fill_per_group(
                fail_histograms["psd"],
                data_m1_psd_fail.geds.energy,
                det_m1_psd_fail,
            )

            data_m2 = data[(data.geds.multiplicity == 2) & good_channel_mask]
            assert mul2_hist is not None
            mul2_hist.fill(
                ak.min(data_m2.geds.energy, axis=-1),
                ak.max(data_m2.geds.energy, axis=-1),
            )

        elif has_spms_coinc:
            pass

    log.info("... convert histograms to lgdo")

    _descriptions = {
        "hit": "all individual HPGe energy deposits in ON channels",
        "mul": "multiplicity-1 events: exactly one ON detector fired",
        "mul_lar": "multiplicity-1 + LAr anti-coincidence",
        "mul_psd": "multiplicity-1 + PSD single-site",
        "mul_lar_psd": "multiplicity-1 + LAr anti-coincidence + PSD single-site",
        "mul2": "multiplicity-2 events: (E_low, E_high) for pairs of ON detectors",
        "fail/lar": "multiplicity-1 events failing the LAr veto",
        "fail/psd": "multiplicity-1 events with valid PSD failing the PSD cut",
    }

    output_dict: dict[str, Struct | Histogram] = {
        cut: Struct(
            {
                g: Histogram(h, attrs={"description": _descriptions[cut]})
                for g, h in per_group.items()
            }
        )
        for cut, per_group in histograms.items()
    }
    if fail_histograms:
        output_dict["fail"] = Struct(
            {
                cut: Struct(
                    {
                        g: Histogram(
                            h, attrs={"description": _descriptions[f"fail/{cut}"]}
                        )
                        for g, h in per_group.items()
                    }
                )
                for cut, per_group in fail_histograms.items()
            }
        )
    if mul2_hist is not None:
        output_dict["mul2"] = Histogram(
            mul2_hist, attrs={"description": _descriptions["mul2"]}
        )
    output = Struct(output_dict)
    lh5.write(output, "pdf", pdf_file, wo_mode="w")
    lh5.write(
        Scalar(int(number_of_primaries)), "nr_sim_events", pdf_file, wo_mode="append"
    )


if __name__ == "__main__":
    main()
