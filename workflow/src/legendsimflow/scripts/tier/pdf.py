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

import awkward as ak
import hist
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
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
    histograms = {
        "hit": h1(),
        "mul": h1(),
        "mul_lar": h1(),
        "mul_psd": h1(),
        "mul_lar_psd": h1(),
        # 2-D: (E_min, E_max) for multiplicity-2 events
        "mul2": hist.new.Reg(6000, 0, 6000).Reg(6000, 0, 6000).Double(),
    }
    fail_histograms = {
        "lar": h1(),
        "psd": h1(),
    }
    log.info("... beginning iteration over cvt file")

    for chunk in iterator:
        data = chunk.view_as("ak")

        # all-good-channel mask: every detector in the event must be usable
        good_channel_mask = ak.all(data.geds.is_good_channel, axis=-1)

        # hit: all energy deposits in good-channel events, no multiplicity requirement
        data_hit = data[good_channel_mask]
        histograms["hit"].fill(ak.flatten(data_hit.geds.energy))

        # multiplicity cut: keep only events where exactly one ON detector fired
        data_m1 = data[(data.geds.multiplicity == 1) & good_channel_mask]
        histograms["mul"].fill(ak.flatten(data_m1.geds.energy))

        # LAr cut: coincident.spms is True when SiPMs detected scintillation light
        # in liquid argon — such events are vetoed
        data_m1_lar = data_m1[~data_m1.coincident.spms]
        histograms["mul_lar"].fill(ak.flatten(data_m1_lar.geds.energy))
        fail_histograms["lar"].fill(
            ak.flatten(data_m1[data_m1.coincident.spms].geds.energy)
        )

        # PSD cut: require valid PSD in data (is_good), simulated A/E observable
        # (has_aoe), and single-site topology (is_single_site). Events where PSD is
        # not valid or not simulated are classified as background and cut.
        def _psd_mask(d):
            return ak.all(
                d.geds.psd.is_good & d.geds.psd.has_aoe & d.geds.psd.is_single_site,
                axis=-1,
            )

        histograms["mul_psd"].fill(ak.flatten(data_m1[_psd_mask(data_m1)].geds.energy))
        histograms["mul_lar_psd"].fill(
            ak.flatten(data_m1_lar[_psd_mask(data_m1_lar)].geds.energy)
        )

        # fail/psd: m1 events with valid and simulated PSD but failing single-site
        psd_fail_mask = ak.all(
            data_m1.geds.psd.is_good & data_m1.geds.psd.has_aoe, axis=-1
        ) & ~_psd_mask(data_m1)
        fail_histograms["psd"].fill(ak.flatten(data_m1[psd_fail_mask].geds.energy))

        # m2: events with exactly two detectors fired — fill 2D histogram with (E_low, E_high)
        data_m2 = data[(data.geds.multiplicity == 2) & good_channel_mask]
        histograms["mul2"].fill(
            ak.min(data_m2.geds.energy, axis=-1), ak.max(data_m2.geds.energy, axis=-1)
        )

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
    output = Struct(
        {
            name: Histogram(h, attrs={"description": _descriptions[name]})
            for name, h in histograms.items()
        }
        | {
            "fail": Struct(
                {
                    name: Histogram(
                        h, attrs={"description": _descriptions[f"fail/{name}"]}
                    )
                    for name, h in fail_histograms.items()
                }
            )
        }
    )
    lh5.write(output, "pdf", pdf_file, wo_mode="w")
    lh5.write(
        Scalar(int(number_of_primaries)), "nr_sim_events", pdf_file, wo_mode="append"
    )


if __name__ == "__main__":
    main()
