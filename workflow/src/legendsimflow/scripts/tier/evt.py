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

from legendsimflow import nersc, patterns
from legendsimflow import reboost as reboost_utils
from legendsimflow.awkward import ak_isin
from legendsimflow.metadata import encode_usability

FORWARD_FIELD_LIST = ["aoe", "drift_time_heuristic", "t0"]
ENERGY_THR_KEV = 25
BUFFER_LEN = "500*MB"
OFF = encode_usability("off")
ON = encode_usability("on")

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

wildcards = args.wildcards
stp_file = patterns.output_simjob_filename(
    args.config, tier="stp", simid=wildcards.simid, jobid=wildcards.jobid
)
hit_file = {
    "opt": args.input.opt_file,
    "hit": args.input.hit_file,
}
evt_file = args.output[0]
log_file = args.log[0]
metadata = args.config.metadata

# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)

log.info("building hit+opt unified TCM")
reboost_utils.build_tcm(hit_file.values(), evt_file)

# test that the evt tcm has the same amount of rows as the stp tcm
assert lh5.read_n_rows("tcm", stp_file) == lh5.read_n_rows("tcm", evt_file)

# get the mapping of detector name to uid
# NOTE: we check on disk because we are not sure which tables were processed in
# the hit tiers
det2uid = {}
for tier in ("opt", "hit"):
    det2uid[tier] = {
        name: uid
        for uid, name in get_remage_detector_uids(
            hit_file[tier], lh5_table="hit"
        ).items()
    }
    msg = f"found mapping name -> uid ({tier} tier): {det2uid[tier]}"
    log.debug(msg)


# little helper to simplify the code below
def _read_hits(tcm_ak, tier, field):
    msg = f"loading {field=} data from {tier=} (file {hit_file[tier]})"
    log.debug(msg)

    return read_data_at_channel_as_ak(
        tcm_ak[tier].table_key,
        tcm_ak[tier].row_in_table,
        hit_file[tier],
        field,
        "hit",
        det2uid[tier],
    )


# iterate over the unified tcm
# NOTE: open mode is append because we will write to the same file
it = lh5.LH5Iterator(evt_file, "tcm", buffer_len=BUFFER_LEN, h5py_open_mode="a")

log.info("begin iterating over TCM")
for chunk in it:
    unified_tcm = chunk.view_as("ak")
    out_table = Table(size=len(unified_tcm))

    # HPGe table
    # ----------
    out_table.add_field("geds", Table(size=len(unified_tcm)))

    # split the unified TCM in two, one for each tier. in this way we will be
    # able to read data from each tier
    tcm = {}
    for tier in ("opt", "hit"):
        mask = ak_isin(unified_tcm.table_key, det2uid[tier].values())
        tcm[tier] = unified_tcm[mask]

    # global fields that are constant over the full events
    # let's take them from the hit tier
    for constant_field in ["run", "period", "evtid"]:
        data = _read_hits(tcm, "hit", constant_field)

        # sanity check
        assert len(data) == len(tcm[tier])

        # replace the awkward missing values with NaN for LH5 compatibility
        data = ak.fill_none(ak.firsts(data, axis=-1), np.nan)
        out_table.add_field(constant_field, Array(np.array(data)))

    # first read usability and energy
    usability = _read_hits(tcm, "hit", "usability")
    energy = _read_hits(tcm, "hit", "energy")

    # we want to only store hits from events in ON and AC detectors and above
    # our energy threshold
    hitsel = (usability != OFF) & (energy > ENERGY_THR_KEV)

    out_table.add_field(
        "geds/is_good_channel", VectorOfVectors(usability[hitsel] == ON)
    )
    out_table.add_field("geds/energy", VectorOfVectors(energy[hitsel]))

    # fields to identify detectors and lookup stuff in the lower tiers
    out_table.add_field("geds/rawid", VectorOfVectors(tcm["hit"].table_key[hitsel]))
    out_table.add_field(
        "geds/hit_idx", VectorOfVectors(tcm["hit"].row_in_table[hitsel])
    )

    # simply forward some fields
    for field in FORWARD_FIELD_LIST:
        field_data = _read_hits(tcm, "hit", field)
        out_table.add_field(f"geds/{field}", VectorOfVectors(field_data[hitsel]))

    # compute multiplicity
    multiplicity = ak.sum(hitsel, axis=-1)
    out_table.add_field("geds/multiplicity", Array(multiplicity))

    # now write down
    lh5.write(out_table, "evt", evt_file, wo_mode="append")
