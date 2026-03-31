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


import awkward as ak
import hist
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
from lgdo import Histogram, Struct, lh5

from legendsimflow import nersc

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

cvt_file = args.input
pdf_file = args.output[0]
log_file = args.log[0]
metadata = args.config.metadata

pdf_file, move2cfs = nersc.make_on_scratch(args.config, cvt_file)

BUFFER_LEN = "500*MB"

# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)

# a single cvt file might be large, need to iterate
iterator = lh5.LH5Iterator(
    str(cvt_file),
    "evt",
    buffer_len=BUFFER_LEN,
)

histograms = {
    "mul_surv": hist.new.Reg(6000, 0, 6000).Double(),
    "hit": hist.new.Reg(6000, 0, 6000).Double(),
    "mul2": hist.new.Reg(6000, 0, 6000).new.Reg(6000, 0, 6000).Double(),
    "mul_lar_surv": hist.new.Reg(6000, 0, 6000).Double(),
}
log.info("... beginning iteration over cvt file")

for chunk in iterator:
    data = chunk.view_as("ak")

    # get hit data (no m1 cut)
    data_hit = data[(ak.all(data.geds.is_good, axis=-1))]

    histograms["hit"].fill(data_hit.geds.energy)

    # get m1 data
    data_m1 = data[(data.geds.multiplicity == 1) & (ak.all(data.geds.is_good, axis=-1))]

    histograms["mul_surv"].fill(data_m1.geds.energy)

    # get m1 data with lar cut
    data_m1_lar = data_m1[(~data_m1.coincident.spms)]
    histograms["mul_lar_surv"].fill(data_m1_lar.geds.energy)

    # get m2 data
    data_m2 = data[(data.geds.multiplicity == 2) & (ak.all(data.geds.is_good, axis=-1))]

    histograms["mul2"].fill(data_m2.geds.energy, data_m2.geds.energy2)

log.info("... convert histograms to lgdo")
histograms_lgdo = Struct({name: Histogram(h) for name, h in histograms.items()})

lh5.write(pdf_file, histograms_lgdo, mode="w")

move2cfs()
