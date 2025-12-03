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
from __future__ import annotations

import logging
from pathlib import Path

import awkward as ak
import h5py
import lgdo
import numpy as np
import pyg4ometry
import pygeomtools
import reboost.hpge.utils
from lgdo import LGDO, lh5

from . import patterns
from .utils import SimflowConfig

log = logging.getLogger(__name__)


def make_output_chunk(chunk: LGDO) -> lgdo.Table:
    out = lgdo.Table(size=len(chunk))

    if "t0" in chunk and isinstance(chunk.t0, lgdo.Array):
        t0 = chunk.t0
    else:
        t0 = lgdo.Array(
            ak.fill_none(ak.firsts(chunk.time.view_as("ak"), axis=-1), 0),
            attrs={"units": "ns"},
        )

    if isinstance(chunk.evtid, lgdo.Array):
        evtid = chunk.evtid
    else:
        evtid = lgdo.Array(
            ak.fill_none(ak.firsts(chunk.evtid.view_as("ak"), axis=-1), 0)
        )

    out.add_field("t0", t0)
    out.add_field("evtid", evtid)

    return out


def write_chunk(
    chunk: lgdo.Table, objname: str, outfile: str, objuid: int, runid: str | None = None
) -> None:
    if not Path(outfile).is_file():
        # create the output file
        lh5.write(lgdo.Struct(), "hit", outfile, wo_mode="write_safe")

    wo_mode = (
        "append_column"
        if objname.strip("/") not in lh5.ls(outfile, "hit/")
        else "append"
    )

    # add runid column to chunk
    # NOTE: is this ok?
    if runid is not None:
        dtype = h5py.string_dtype(encoding="utf-8", length=len(runid))
        runid = np.full(len(chunk), fill_value=runid, dtype=dtype)
        chunk.add_field("runid", lgdo.Array(runid))

    lh5.write(
        chunk,
        objname,
        outfile,
        wo_mode=wo_mode,
    )

    if wo_mode == "append_column":
        if "hit/__by_uid__" not in lh5.ls(outfile, "hit/"):
            log.debug("creating hit/__by_uid__ folder")
            lh5.write(lgdo.Struct(), "hit/__by_uid__", outfile)

        msg = f"creating soft link hit/__by_uid__/det{objuid} -> {objname}"
        log.debug(msg)
        with h5py.File(outfile, "r+") as f:
            # create uid -> det_name symlink
            f[f"hit/__by_uid__/det{objuid}"] = h5py.SoftLink(objname)
            # updated the struct datatype attribute by adding the new symlink
            dt = f["hit/__by_uid__"].attrs.pop("datatype")
            fields = [*lgdo.lh5.datatype.get_struct_fields(dt), f"det{objuid}"]
            f["hit/__by_uid__"].attrs["datatype"] = (
                "struct{" + ",".join(sorted(fields)) + "}"
            )


def get_sensvols(geom, det_type: str | None = None) -> list[str]:
    sensvols = pygeomtools.detectors.get_all_sensvols(geom)
    if det_type is not None:
        return [k for k, v in sensvols.items() if v.detector_type == det_type]
    return list(sensvols.keys())


def load_hpge_dtmaps(
    config: SimflowConfig, det_name: str, runid: str
) -> dict[str, reboost.hpge.utils.HPGeScalarRZField]:
    hpge_dtmap_file = patterns.output_dtmap_merged_filename(
        config,
        runid=runid,
    )

    if len(lh5.ls(hpge_dtmap_file, f"{det_name}/drift_time_*")) >= 2:
        log.debug("loading drift time maps")
        dt_map = {}
        for angle in ("000", "045"):
            dt_map[angle] = reboost.hpge.utils.get_hpge_scalar_rz_field(
                hpge_dtmap_file, det_name, f"drift_time_{angle}_deg"
            )
    else:
        msg = (
            f"no valid time maps found for {det_name} in {hpge_dtmap_file}, "
            "drift time will be set to NaN"
        )
        log.warning(msg)
        dt_map = None


def get_remage_hit_range(
    stp_file: str | Path, det_name: str, uid: int, evt_idx_range: list[int]
) -> tuple[int]:
    # load TCM, to be used to chunk the event statistics according to the run partitioning
    tcm = lh5.read_as("tcm", stp_file, library="ak")

    # ask the TCM which rows we should read from the hit table
    tcm_part = tcm[evt_idx_range[0] : evt_idx_range[1]]
    entry_list = ak.flatten(tcm_part[tcm_part.table_key == uid].row_in_table).to_list()

    if len(entry_list) > 0:
        assert list(range(entry_list[0], entry_list[-1] + 1)) == entry_list

        msg = (
            f"hits with indices in [{entry_list[0]}, {entry_list[-1]}] "
            "recorded in the events belonging to this partition"
        )
        log.debug(msg)

        i_start = entry_list[0]
        n_entries = entry_list[-1] - entry_list[0]

    else:
        msg = (
            f"no hits recorded in {det_name} in the events belonging to this partition"
        )
        log.warning(msg)

        i_start = 0
        n_entries = None

    return i_start, n_entries


def hpge_corrected_dt_heuristic(
    chunk: ak.Array,
    dt_map: reboost.hpge.utils.HPGeScalarRZField,
    det_loc: pyg4ometry.gdml.Defines.Position,
) -> ak.Array:
    _phi = np.arctan2(
        chunk.yloc * 1000 - det_loc.eval()[1],
        chunk.xloc * 1000 - det_loc.eval()[0],
    )

    _drift_time = {}
    for angle, _map in dt_map.items():
        _drift_time[angle] = reboost.hpge.psd.drift_time(
            chunk.xloc,
            chunk.yloc,
            chunk.zloc,
            _map,
            coord_offset=det_loc,
        ).view_as("ak")

    _drift_time_corr = (
        _drift_time["045"]
        + (_drift_time["000"] - _drift_time["045"]) * (1 - np.cos(4 * _phi)) / 2
    )

    return reboost.hpge.psd.drift_time_heuristic(_drift_time_corr, chunk.edep)
