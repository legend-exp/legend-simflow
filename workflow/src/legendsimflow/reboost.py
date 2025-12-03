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
import pygeomtools
from lgdo import LGDO, lh5

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


def write_chunk(chunk, objname, outfile, objuid, runid):
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
