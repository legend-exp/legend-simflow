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
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import numpy as np
from lgdo import Array, Table, VectorOfVectors, lh5
from reboost.core import read_data_at_channel_as_ak
from reboost.utils import get_remage_detector_uids

from legendsimflow import nersc
from legendsimflow import reboost as reboost_utils
from legendsimflow.metadata import encode_usability

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

stp_file = args.input.stp_file
opt_file = args.input.opt_file
hit_file = args.input.hit_file
evt_file = args.output[0]
log_file = args.log[0]
metadata = args.config.metadata
buffer_len = args.params.buffer_len


# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)

log.info("building hit+opt unified TCM")
reboost_utils.build_tcm([hit_file, opt_file], evt_file)

# test that the evt tcm has the same amount of rows as the stp tcm
assert lh5.read_n_rows("tcm", stp_file) == lh5.read_n_rows("tcm", evt_file)

log.info("extracting remage uid map")

# get the mapping of detector name to uid
log.info(type(hit_file))
mapping = {
    name: uid
    for uid, name in get_remage_detector_uids(hit_file, lh5_table="hit").items()
}

# fields to forward
field_list = ["aoe", "drift_time_heuristic", "energy", "t0"]

log.info("begin iterating over unified TCM")

off = encode_usability("off")
on = encode_usability("ac")

# iterate
it = lh5.LH5Iterator(hit_file, "tcm", buffer_len=buffer_len)

for tcm_chunk in it:
    # convert tcm to ak
    tcm_tmp = tcm_chunk.view_as("ak")

    # create output
    out = Table(size=len(tcm_tmp))

    # first read usability
    usability = read_data_at_channel_as_ak(
        tcm_tmp.table_key, tcm_tmp.row_in_table, hit_file, "usability", "hit", mapping
    )

    out.add_field("usability", VectorOfVectors(usability[usability != off]))

    # constant over the full events
    for constant_field in ["run", "period", "evtid"]:
        data = read_data_at_channel_as_ak(
            tcm_tmp.table_key,
            tcm_tmp.row_in_table,
            hit_file,
            constant_field,
            "hit",
            mapping,
        )

        if not ak.all(data == data[:, 0]):
            msg = f"All values of {constant_field} in the event must be the same"
            raise ValueError(msg)

        data = ak.fill_none(ak.firsts(data, axis=-1), np.nan)
        out.add_field(constant_field, Array(np.array(data)))

    # now lets add all channel dependent fields
    usability_no_off = usability[usability != off]

    out.add_field("channel", VectorOfVectors(tcm_tmp.table_key[usability != off]))
    out.add_field("row", VectorOfVectors(tcm_tmp.row_in_table[usability != off]))

    # add other fields
    for field in field_list:
        field_data = read_data_at_channel_as_ak(
            tcm_tmp.table_key, tcm_tmp.row_in_table, hit_file, field, "hit", mapping
        )

        # filter based on usability
        field_data = field_data[usability != off]

        out.add_field(field, VectorOfVectors(field_data))

    # compute multiplicity and check if the event is good
    energy = out.energy.view_as("ak")

    is_good = ak.all(usability_no_off[energy > 25] == on, axis=-1)

    multiplicity = ak.sum(energy > 25, axis=-1)
    multiplicity_good = ak.sum(energy[usability_no_off == 0] > 25, axis=-1)

    out.add_field("multiplicity", Array(multiplicity))
    out.add_field("is_good", Array(is_good))
    out.add_field("multiplicity_good", Array(multiplicity_good))

    # now write
    lh5.write(out, "evt", evt_file, wo_mode="append")
