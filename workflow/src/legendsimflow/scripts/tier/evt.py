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

FORWARD_FIELD_LIST = ["aoe", "drift_time_heuristic", "t0"]
ENERGY_THR_KEV = 25
OFF = encode_usability("off")
ON = encode_usability("on")
BUFFER_LEN = "500*MB"

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

stp_file = args.input.stp_file
opt_file = args.input.opt_file
hit_file = args.input.hit_file
evt_file = args.output[0]
log_file = args.log[0]
metadata = args.config.metadata

# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)

log.info("building hit+opt unified TCM")
reboost_utils.build_tcm([hit_file, opt_file], evt_file)

# test that the evt tcm has the same amount of rows as the stp tcm
assert lh5.read_n_rows("tcm", stp_file) == lh5.read_n_rows("tcm", evt_file)

# get the mapping of detector name to uid
detname2uid = {
    name: uid
    for uid, name in get_remage_detector_uids(hit_file, lh5_table="hit").items()
}


def _read_hits(tcm_ak, field):
    return read_data_at_channel_as_ak(
        tcm_ak.table_key, tcm_ak.row_in_table, hit_file, field, "hit", detname2uid
    )


# iterate
it = lh5.LH5Iterator(hit_file, "tcm", buffer_len=BUFFER_LEN)

log.info("begin iterating over TCM")
for tcm_chunk in it:
    tcm_ak = tcm_chunk.view_as("ak")
    out_table = Table(size=len(tcm_ak))

    # global fields that are constant over the full events
    for constant_field in ["run", "period", "evtid"]:
        data = _read_hits(tcm_ak, constant_field)

        # sanity check
        assert ak.all(data == data[:, 0])

        # replace the awkward missing values with NaN for LH5 compatibility
        data = ak.fill_none(ak.firsts(data, axis=-1), np.nan)
        out_table.add_field(constant_field, Array(np.array(data)))

    # HPGe table
    # ----------
    out_table.add_field("geds", Table(size=len(tcm_ak)))

    # first read usability and energy
    usability = _read_hits(tcm_ak, "usability")
    energy = _read_hits(tcm_ak, "energy")

    # we want to only store hits from events in ON and AC detectors and above
    # our energy threshold
    is_good_hit = (usability != OFF) & (energy > ENERGY_THR_KEV)

    out_table.add_field("geds/usability", VectorOfVectors(usability[is_good_hit]))
    out_table.add_field("geds/energy", VectorOfVectors(energy[is_good_hit]))

    # fields to identify detectors and lookup stuff in the lower tiers
    out_table.add_field("geds/rawid", VectorOfVectors(tcm_ak.table_key[is_good_hit]))
    out_table.add_field(
        "geds/hit_idx", VectorOfVectors(tcm_ak.row_in_table[is_good_hit])
    )

    # simply forward some fields
    for field in FORWARD_FIELD_LIST:
        field_data = _read_hits(tcm_ak, field)
        out_table.add_field(f"geds/{field}", VectorOfVectors(field_data[is_good_hit]))

    # compute multiplicity
    multiplicity = ak.sum(is_good_hit, axis=-1)
    out_table.add_field("geds/multiplicity", Array(multiplicity))

    # now write down
    lh5.write(out_table, "evt", evt_file, wo_mode="append")
