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

    cvt_file = nersc.dvs_ro(config, args.cvt_file)
    pdf_file, move2cfs = nersc.make_on_scratch(config, args.pdf_file)
    log_file = args.log_file
    simid = args.simid
    experiment = config.experiment

    BUFFER_LEN = "500*MB"

    log = ldfs.utils.build_log(config.metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "tier-pdf", parser, args)

    msg = f"... reading data from {cvt_file}"
    log.info(msg)

    # a single cvt file might be large, need to iterate
    iterator = lh5.LH5Iterator(
        str(cvt_file),
        "evt",
        buffer_len=BUFFER_LEN,
    )

    simconfig_block = config.metadata.simprod.config.tier.stp[experiment].simconfig[
        simid
    ]
    number_of_primaries = (
        simconfig_block.primaries_per_job * simconfig_block.number_of_jobs
    )

    msg = f"... extracted number of primaries from simconfig {number_of_primaries}"
    log.info(msg)

    # 1 keV/bin over 0-6000 keV; boost-histogram only provides Double (float64) storage
    histograms = {
        "mul_surv": hist.new.Reg(6000, 0, 6000).Double(),
        "hit": hist.new.Reg(6000, 0, 6000).Double(),
        "mul2": hist.new.Reg(6000, 0, 6000).Reg(6000, 0, 6000).Double(),
        "mul_lar_surv": hist.new.Reg(6000, 0, 6000).Double(),
    }
    log.info("... beginning iteration over cvt file")

    for chunk in iterator:
        data = chunk.view_as("ak")

        # get hit data (no m1 cut)
        data_hit = data[(ak.all(data.geds.is_good_channel, axis=-1))]

        histograms["hit"].fill(ak.flatten(data_hit.geds.energy))

        # get m1 data
        data_m1 = data[
            (data.geds.multiplicity == 1) & (ak.all(data.geds.is_good_channel, axis=-1))
        ]

        histograms["mul_surv"].fill(ak.flatten(data_m1.geds.energy))

        # get m1 data with lar cut
        data_m1_lar = data_m1[(~data_m1.coincident.spms)]
        histograms["mul_lar_surv"].fill(ak.flatten(data_m1_lar.geds.energy))

        # get m2 data
        data_m2 = data[
            (data.geds.multiplicity == 2) & (ak.all(data.geds.is_good_channel, axis=-1))
        ]

        histograms["mul2"].fill(
            ak.min(data_m2.geds.energy, axis=-1), ak.max(data_m2.geds.energy, axis=-1)
        )

    log.info("... convert histograms to lgdo")
    histograms_lgdo = Struct({name: Histogram(h) for name, h in histograms.items()})

    output = Struct({"pdfs": histograms_lgdo, "n_prim": Scalar(number_of_primaries)})
    lh5.write(output, "pdf", pdf_file, wo_mode="w")

    move2cfs()


if __name__ == "__main__":
    main()
