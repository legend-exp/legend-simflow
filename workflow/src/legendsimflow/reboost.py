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
from collections.abc import Iterable, Mapping
from pathlib import Path

import awkward as ak
import h5py
import lgdo
import numpy as np
import pyg4ometry
import pygama.evt
import pygeomtools
import reboost.hpge.utils
from lgdo import LGDO, lh5

from . import patterns
from .utils import SimflowConfig

log = logging.getLogger(__name__)


def make_output_chunk(chunk: LGDO) -> lgdo.Table:
    """Prepare output detector table chunk for the hit tier.

    Note
    ----
    This function will be moved to :mod:`reboost`.
    """
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


def write_chunk(chunk: lgdo.Table, objname: str, outfile: str, objuid: int) -> None:
    """Write detector table chunks for the hit tier to disk.

    Note
    ----
    This function will be moved to :mod:`reboost`.
    """
    if not Path(outfile).is_file():
        # create the output file
        lh5.write(lgdo.Struct(), "hit", outfile, wo_mode="write_safe")

    wo_mode = (
        "append_column"
        if objname.strip("/") not in lh5.ls(outfile, "hit/")
        else "append"
    )

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

        link_name = f"det{objuid:03}"

        msg = f"creating soft link hit/__by_uid__/{link_name} -> {objname}"
        log.debug(msg)
        with h5py.File(outfile, "r+", locking=False) as f:
            # create uid -> det_name symlink
            f[f"hit/__by_uid__/{link_name}"] = h5py.SoftLink(objname)
            # updated the struct datatype attribute by adding the new symlink
            dt = f["hit/__by_uid__"].attrs.pop("datatype")
            fields = [*lgdo.lh5.datatype.get_struct_fields(dt), link_name]
            f["hit/__by_uid__"].attrs["datatype"] = (
                "struct{" + ",".join(sorted(fields)) + "}"
            )


def get_senstables(
    geom: pyg4ometry.geant4.Registry, det_type: str | None = None
) -> list[str]:
    sensvols = pygeomtools.detectors.get_all_senstables(geom)
    if det_type is not None:
        return [k for k, v in sensvols.items() if v.detector_type == det_type]
    return list(sensvols.keys())


def load_hpge_dtmaps(
    config: SimflowConfig, det_name: str, runid: str
) -> dict[str, reboost.hpge.utils.HPGeRZField]:
    """Load HPGe drift time maps from disk.

    Automatically finds and loads drift time maps for crystal axes <100> <110>.

    Note
    ----
    This function will be moved to :mod:`reboost`.
    """
    hpge_dtmap_file = patterns.output_dtmap_merged_filename(
        config,
        runid=runid,
    )

    if len(lh5.ls(hpge_dtmap_file, f"{det_name}/drift_time_*")) >= 2:
        log.debug("loading drift time maps")
        dt_map = {}
        for angle in ("000", "045"):
            dt_map[angle] = reboost.hpge.utils.get_hpge_rz_field(
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
    tcm: ak.Array, det_name: str, uid: int, evt_idx_range: list[int]
) -> tuple[int]:
    """Extract the range of remage output rows for an event range.

    Queries the remage TCM (stored below ``/tcm`` in `stp_file`) with the input
    `evt_idx_range = [i, j]` to extract the first and last index of rows (hits)
    in the `det_name` detector table that correspond to the input event range.
    Returns the start index and number of rows to read after it as a tuple.

    Parameters
    ----------
    tcm
        Time-coincidence map.
    det_name
        name of the detector table in `stp_file`.
    uid
        remage unique identifier for detector `det_name`.
    evt_idx_range
        `[first, last]` (i.e. `first` included, `last` included) index of
        events of interest present in the remage output file. Only positive
        indices are supported.
    """

    if (evt_idx_range[0] < 0) or (evt_idx_range[0] < 0):
        msg = "Only positive indices are supported"
        raise ValueError(msg)

    # add one for inclusive slicing
    tcm_part = tcm[evt_idx_range[0] : evt_idx_range[1] + 1]

    entry_list = ak.flatten(tcm_part[tcm_part.table_key == uid].row_in_table).to_list()

    assert entry_list == sorted(entry_list)

    if len(entry_list) > 0:
        assert list(range(entry_list[0], entry_list[-1] + 1)) == entry_list

        msg = (
            f"hits with indices in [{entry_list[0]}, {entry_list[-1]}] "
            "recorded in the events belonging to this partition"
        )
        log.debug(msg)

        i_start = entry_list[0]
        n_entries = entry_list[-1] - entry_list[0] + 1

    else:
        msg = (
            f"no hits recorded in {det_name} in the events belonging to this partition"
        )
        log.warning(msg)

        i_start = 0
        n_entries = None

    return i_start, n_entries


def hpge_corrected_drift_time(
    chunk: ak.Array,
    dt_map: reboost.hpge.utils.HPGeRZField,
    det_loc: pyg4ometry.gdml.Defines.Position,
) -> ak.Array:
    """HPGe drift time heuristic corrected for crystal axis effects.

    Note
    ----
    This function will be moved to :mod:`reboost`.
    """
    phi = np.arctan2(
        chunk.yloc * 1000 - det_loc.eval()[1],
        chunk.xloc * 1000 - det_loc.eval()[0],
    )

    drift_time = {}
    for angle, _map in dt_map.items():
        drift_time[angle] = reboost.hpge.psd.drift_time(
            chunk.xloc,
            chunk.yloc,
            chunk.zloc,
            _map,
            coord_offset=det_loc,
        ).view_as("ak")

    return (
        drift_time["045"]
        + (drift_time["000"] - drift_time["045"]) * (1 - np.cos(4 * phi)) / 2
    )


def hpge_max_current_cal(
    edep: ak.Array,
    drift_time: ak.Array,
    currmod_pars: Mapping,
) -> ak.Array:
    """Calculate the maximum of the current pulse.

    Parameters
    ----------
    edep
        energy deposited at each step.
    drift_time
        drift time of each energy deposit.
    currmod_pars
        dictionary storing the parameters of the current model (see
        :func:`reboost.hpge.psd.get_current_template`)
    """
    # current pulse template domain in ns (step is 1 ns)
    t_domain = {"low": -1000, "high": 4000, "step": 1}

    # instantiate the template
    a_tmpl = reboost.hpge.psd.get_current_template(
        **t_domain,
        mean_aoe=1,  # set the maximum of the template to unity, so the A/E will be calibrated
        **currmod_pars,
    )
    # and calculate the maximum current
    return reboost.hpge.psd.maximum_current(
        edep,
        drift_time,
        template=a_tmpl,
        times=np.arange(t_domain["low"], t_domain["high"]),
    )


def build_tcm(
    hit_files: str | Path | Iterable[str | Path], out_file: str | Path
) -> None:
    """Re-create the TCM table from remage.

    Use remage fields `evtid` and `t0` (the latter is assumed to be in
    nanoseconds) to build coincidences. The settings are identical to the
    remage built-in TCM settings.

    Note
    ----
    This function will be moved to :mod:`reboost`.
    """
    if isinstance(hit_files, str | Path):
        hit_files = [hit_files]

    # use tables keyed by UID in the __by_uid__ group.  in this way, the
    # TCM will index tables by UID.  the coincidence criterium is based
    # on Geant4 event identifier and time of the hits
    # NOTE: uses the same time window as in build_hit() reshaping
    pygama.evt.build_tcm(
        [(f, r"hit/__by_uid__/*") for f in hit_files],  # input_tables
        ["evtid", "t0"],  # coin_cols
        hash_func=r"(?<=hit/__by_uid__/det)\d+",
        coin_windows=[0, 10_000],  # in ns
        out_file=out_file,
        wo_mode="write_safe",
    )
