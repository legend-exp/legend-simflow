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
import re
from collections.abc import Mapping
from pathlib import Path

import awkward as ak
import h5py
import lgdo
import numba as nb
import numpy as np
import pyg4ometry
import pygeomtools
import reboost.hpge.utils
import reboost.units
from lgdo import LGDO, lh5

from . import patterns
from .utils import SimflowConfig

log = logging.getLogger(__name__)


def get_remage_detector_uids(h5file: str | Path, *, lh5_table: str = "stp") -> dict:
    """Get mapping of detector names to UIDs from a remage output file.

    The remage LH5 output files contain a link structure that lets the user
    access detector tables by UID. For example:

    .. code-block:: text

        ├── stp · struct{det1,det2,optdet1,optdet2,scint1,scint2}
        └── __by_uid__ · struct{det001,det002,det011,det012,det101,det102}
            ├── det001 -> /stp/scint1
            ├── det002 -> /stp/scint2
            ├── det011 -> /stp/det1
            ├── det012 -> /stp/det2
            ├── det101 -> /stp/optdet1
            └── det102 -> /stp/optdet2

    This function analyzes this structure and returns:

    .. code-block:: text

        {1: 'scint1',
         2: 'scint2',
         11: 'det1',
         12: 'det2',
         101: 'optdet1',
         102: 'optdet2'g

    Parameters
    ----------
    h5file
        path to remage output file.
    """
    if isinstance(h5file, Path):
        h5file = h5file.as_posix()

    out = {}
    with h5py.File(h5file, "r", locking=False) as f:
        g = f[f"/{lh5_table}/__by_uid__"]
        # loop over links
        for key in g:
            # is this a link?
            link = g.get(key, getlink=True)
            if isinstance(link, h5py.SoftLink):
                m = re.fullmatch(r"det(\d+)", key)
                if m is None:
                    msg = rf"'{key}' is not formatted as expected, i.e. 'det(\d+)', skipping"
                    log.warning(msg)
                    continue

                # get the name of the link target without trailing groups (to
                # i.e. remove /stp)
                name = link.path.split("/")[-1]

                out[int(m.group(1))] = name
    return out


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
) -> dict[str, reboost.hpge.utils.HPGeRZField] | None:
    """Load HPGe drift time maps from disk.

    Automatically finds and loads drift time maps for crystal axes <100> <110>.
    If no map is found, None is returned.

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
                hpge_dtmap_file,
                det_name,
                f"drift_time_{angle}_deg",
                bounds_error=False,
            )
    else:
        msg = (
            f"no valid time maps found for {det_name} in {hpge_dtmap_file}, "
            "drift time will be set to NaN"
        )
        log.warning(msg)
        dt_map = None

    return dt_map


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
        n_entries = 0

    return i_start, n_entries


def hpge_corrected_drift_time(
    chunk: ak.Array,
    dt_map: dict[str, reboost.hpge.utils.HPGeRZField],
    det_loc: pyg4ometry.gdml.Defines.Position,
) -> ak.Array:
    """HPGe drift time heuristic corrected for crystal axis effects.

    Note
    ----
    This function will be moved to :mod:`reboost`.
    """
    # Convert det_loc to pint Quantity
    det_loc_pint = reboost.units.pg4_to_pint(det_loc)

    # Use reboost.units to get conversion factors for chunk coordinates
    # This handles the case when chunk has units attached (with_units=True)
    xloc_conv = reboost.units.units_convfact(chunk.xloc, det_loc_pint.units)
    yloc_conv = reboost.units.units_convfact(chunk.yloc, det_loc_pint.units)

    # Unwrap LGDO/pint if present
    xloc, _ = reboost.units.unwrap_lgdo(chunk.xloc)
    yloc, _ = reboost.units.unwrap_lgdo(chunk.yloc)

    phi = np.arctan2(
        yloc * yloc_conv - det_loc_pint[1].m,
        xloc * xloc_conv - det_loc_pint[0].m,
    )

    drift_time = {}
    for angle, _map in dt_map.items():
        drift_time[angle] = reboost.hpge.psd.drift_time(
            chunk.xloc,
            chunk.yloc,
            chunk.zloc,
            _map,
            coord_offset=det_loc,
        )

    return (
        drift_time["045"]
        + (drift_time["000"] - drift_time["045"]) * (1 - np.cos(4 * phi)) / 2
    )


def hpge_max_current(
    edep: ak.Array,
    drift_time: ak.Array,
    currmod_pars: Mapping,
    **kwargs,
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
    kwargs
        forwarded to :func:`reboost.hpge.psd.maximum_current`.
    """
    # current pulse template domain in ns (step is 1 ns)
    t_domain = {"low": -1000, "high": 4000, "step": 1}

    # instantiate the template
    a_tmpl, times = reboost.hpge.psd.get_current_template(
        **t_domain,
        mean_aoe=1,  # set the maximum of the template to unity, so the A/E will be calibrated
        **currmod_pars,
    )
    # and calculate the maximum current
    return reboost.hpge.psd.maximum_current(
        edep,
        drift_time,
        template=a_tmpl,
        times=times,
        **kwargs,
    )


@nb.njit(cache=True)
def _cluster_photoelectrons_flat(
    offsets: np.ndarray,
    t: np.ndarray,
    a: np.ndarray,
    thr: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Numba-accelerated clustering kernel for innermost list level.

    Parameters
    ----------
    offsets
        1D int64 array of list offsets (length = n_lists + 1).
    t
        1D array of sorted times for all elements.
    a
        1D array of amplitudes corresponding to times.
    thr
        Maximum time span within a cluster.

    Returns
    -------
    out_t
        Clustered times (first time in each cluster).
    out_a
        Clustered amplitudes (sum of amplitudes in each cluster).
    counts
        Number of clusters per original list.
    """
    n_lists = offsets.size - 1
    n = t.size

    out_t = np.empty(n, dtype=t.dtype)
    out_a = np.empty(n, dtype=a.dtype)
    counts = np.empty(n_lists, dtype=np.int64)

    out_k = 0
    for ilist in range(n_lists):
        s = offsets[ilist]
        e = offsets[ilist + 1]

        ncl = 0
        i = s
        while i < e:
            t0 = t[i]
            asum = 0.0
            j = i
            while j < e and (t[j] - t0) <= thr:
                asum += a[j]
                j += 1

            out_t[out_k] = t0
            out_a[out_k] = asum
            out_k += 1
            ncl += 1
            i = j

        counts[ilist] = ncl

    return out_t[:out_k], out_a[:out_k], counts


def _listoffset_chain(
    layout: ak.contents.Content,
) -> tuple[list[np.ndarray], ak.contents.NumpyArray]:
    """Extract the chain of offsets from nested ListOffsetArrays.

    Parameters
    ----------
    layout
        An awkward array layout.

    Returns
    -------
    offsets_chain
        List of np.int64 arrays, one per nested list depth,
        from outermost to innermost.
    content_numpy_layout
        The final NumpyArray content node.
    """
    offsets_chain = []
    node = layout
    while isinstance(node, ak.contents.ListOffsetArray):
        offsets_chain.append(np.asarray(node.offsets, dtype=np.int64))
        node = node.content
    if not isinstance(node, ak.contents.NumpyArray):
        msg = (
            "expected ListOffsetArray(s) ending in NumpyArray; got "
            f"{type(node).__name__}"
        )
        raise TypeError(msg)
    return offsets_chain, node


def cluster_photoelectrons(
    times: ak.Array, amps: ak.Array, thr: float
) -> tuple[ak.Array, ak.Array]:
    """Cluster photoelectrons within the instrument time resolution.

    Clusters hits at axis=-1 (innermost lists) such that within each cluster
    the time span (last_time - first_time) does not exceed `thr`. This is
    useful for combining photoelectrons that arrive within the time resolution
    of the detector, treating them as a single detected event.

    The output time is the first time of each cluster; the amplitude is the
    sum of all amplitudes in the cluster.

    Parameters
    ----------
    times
        Awkward array of hit times. Must be sorted in ascending order within
        each innermost list. Sorting is the caller's responsibility; unsorted
        input produces undefined behavior.
    amps
        Awkward array of amplitudes corresponding to times. Must have the same
        structure (nesting depth and list lengths) as `times`.
    thr
        Maximum time span within a cluster (e.g., the detector time resolution).

    Returns
    -------
    clustered_times
        Awkward array with the same nesting structure, containing the first
        time of each cluster.
    clustered_amps
        Awkward array with the same nesting structure, containing the summed
        amplitude of each cluster.

    Raises
    ------
    ValueError
        If `times` and `amps` have different nesting depths or different
        numbers of elements.

    Examples
    --------
    >>> times = ak.Array([[0.0, 0.6, 1.1, 1.4, 2.3]])
    >>> amps = ak.Array([[1.0, 2.0, 3.0, 4.0, 5.0]])
    >>> t_out, a_out = cluster_photoelectrons(times, amps, thr=1.0)
    >>> ak.to_list(t_out)
    [[0.0, 1.1, 2.3]]
    >>> ak.to_list(a_out)
    [[3.0, 7.0, 5.0]]
    """
    times_p = ak.to_packed(times)
    amps_p = ak.to_packed(amps)

    lt = ak.to_layout(times_p)
    la = ak.to_layout(amps_p)

    offsets_chain_t, tnode = _listoffset_chain(lt)
    offsets_chain_a, anode = _listoffset_chain(la)

    # minimal shape sanity: same nesting depth for list structure
    if len(offsets_chain_t) != len(offsets_chain_a):
        msg = "times and amps do not have the same list nesting depth"
        raise ValueError(msg)

    # verify that offsets match at each nesting level
    for i, (off_t, off_a) in enumerate(
        zip(offsets_chain_t, offsets_chain_a, strict=True)
    ):
        if not np.array_equal(off_t, off_a):
            msg = f"times and amps have mismatched list lengths at nesting level {i}"
            raise ValueError(msg)

    # innermost offsets define the axis=-1 lists that must not mix
    inner_offsets = offsets_chain_t[-1]

    # 1D content buffers
    t = np.asarray(tnode.data)
    a = np.asarray(anode.data)

    out_t, out_a, inner_counts = _cluster_photoelectrons_flat(inner_offsets, t, a, thr)

    # rebuild: first make clustered axis=-1 lists
    cur_t = ak.unflatten(out_t, inner_counts)
    cur_a = ak.unflatten(out_a, inner_counts)

    # then rebuild outer list levels (working outward)
    # each parent list has counts = diff(parent_offsets)
    for parent_offsets in reversed(offsets_chain_t[:-1]):
        parent_counts = np.diff(parent_offsets).astype(np.int64)
        cur_t = ak.unflatten(cur_t, parent_counts)
        cur_a = ak.unflatten(cur_a, parent_counts)

    return cur_t, cur_a


def smear_photoelectrons(
    array: ak.Array, fwhm_in_pe: float, rng: np.random.Generator = None
) -> ak.Array:
    """Smear photoelectron pulse amplitudes.

    Returns an array of gaussian distributed single-photoelectron amplitudes
    with the same shape of the input `array`.
    """
    if rng is None:
        rng = np.random.default_rng()

    counts = ak.num(array)
    flat = rng.normal(loc=1, scale=fwhm_in_pe / 2.35482, size=ak.sum(counts))
    flat = np.where(flat < 0, 0, flat)

    return ak.unflatten(flat, counts)


def gauss_smear(arr_true: ak.Array, arr_reso: ak.Array) -> ak.Array:
    """Smear values with expected resolution.

    Samples from gaussian and shifts negative values to a fixed, tiny positive
    value.
    """
    arr_smear = reboost.math.stats.gaussian_sample(
        arr_true,
        arr_reso,
    )

    # energy can't be negative as a result of smearing
    return ak.where((arr_smear <= 0) & (arr_true >= 0), np.finfo(float).tiny, arr_smear)

def get_forced_trigger_library(
    evt_file_path: str | Path,
    win_len_ns: float = 6_000,
    min_sep_ns: float = 6_000,
    win_range_ns: tuple[float, float] = (1_000, 44_000),
) -> ak.Array:
    """Extract a library of forced trigger events to be used
    in correcting the SiPM pe and times with random coincidences.

    This reformats the data to make use of several windows within a waveform to build forced trigger events
    and then stores the number of pe and times for each SiPM channel from that window (to be used for corrections).

    Using this approach, we can maximize the number of forced-trigger-like events,
    that are not actually generated from forced trigger data.

    Parameters
    ----------
    evt_file_path
        Path to the event tier data.
    win_len
        Window length of a "forced trigger", in nanoseconds.
    min_sep
        Minimal separation time between two windows in a trace, in nanoseconds.
    win_range_ns
        Available window range in a waveform in data.
        Depends on which data is used (forced trigger data allows full waveform, while HPGe-/LAr-triggered data only allows first half).

    Returns
    -------
    Array with fields "npe", the number of pe per SiPM and per hit,
    "t0", the time relative to the start of a window in the trace, per SiPM and per hit (makes sure that t0 are between 0 and win_len time),
    and "rawids" the SiPM channel numbers.
    """

    npe = None
    t0 = None
    rawids = None

    # determine window starts and ends
    starts = []
    current_start = win_range_ns[0]
    while current_start + win_len_ns <= win_range_ns[1]:
        starts.append(current_start)
        current_start += win_len_ns + min_sep_ns  # move to next window
    ends = [s + win_len_ns for s in starts]

    files = Path(evt_file_path).glob("*")

    for file in files:
        spms = lh5.read(
            "evt/",
            file,
            field_mask=["spms/energy", "spms/t0", "spms/rawid"],
        ).view_as("ak")

        rawids_tmp = spms.spms.rawid[0]

        if rawids is not None and not ak.all(rawids == rawids_tmp):
            msg = "rawid should be the same in all cases"
            raise ValueError(msg)

        for wstart, wend in zip(starts, ends, strict=True):
            tmsk = (spms.spms.t0 > wstart) & (spms.spms.t0 <= wend)
            npe_tmp = spms.spms.energy[tmsk]

            # times between 0 and 6 µs
            t0_tmp = spms.spms.t0[tmsk] - wstart

            npe = npe_tmp if npe is None else ak.concatenate((npe, npe_tmp))

            t0 = t0_tmp if t0 is None else ak.concatenate((t0, t0_tmp))

        rawids = rawids_tmp

    return ak.Array(
        {
            "npe": npe,
            "t0": t0,
            "rawid": np.vstack([rawids] * len(npe)),
        }
    )


def add_forced_triggers(
    sipm_uid: int, rc_data: ak.Array, times: ak.Array, amps: ak.Array
) -> tuple[ak.Array, ak.Array]:
    """Apply a random coincidence correction to simulated number of pe and times.

    Parameters
    ----------
    sipm_uid
        The unique identifier (uid) for a SiPM channel.
    rc_data
        Random coincidence data library with fields "npe", "t0", "rawid".
    times
        Awkward array of hit times.
    amps
        Awkward array of amplitudes corresponding to times.
    Returns
    -------
    corr_t
        Awkward array containing simulated hit times, and the random coincidence hit times.
    corr_a
        Awkward array containing simulated hit amplitues, and the random coincidence hit amplitudes.
    """

    rc_data_ch = rc_data[rc_data.rawid == sipm_uid]

    # Check that times and amps have same length
    if len(times) != len(amps):
        msg = f"times length ({len(times)}) != amps length ({len(amps)})"
        raise ValueError(msg)

    # Check that rc_data.t0 and rc_data.npe have same length
    if len(rc_data_ch.t0) != len(rc_data_ch.npe):
        msg = f"rc_data.t0 length ({len(rc_data_ch.t0)}) != rc_data.npe length ({len(rc_data_ch.npe)})"
        raise ValueError(msg)

    if len(rc_data_ch) < len(amps):
        msg = f"Random coincidence data length ({len(rc_data_ch)}) is smaller than simulated data length ({len(amps)})"
        raise ValueError(msg)

    if len(rc_data_ch) > len(amps):
        rng = np.random.default_rng()
        indices = rng.choice(len(rc_data_ch), size=len(amps), replace=False)
        rc_data_ch = rc_data_ch[
            np.sort(indices)
        ]  # sorting to keep temporal order of forced trigger data, could be removed if unwanted.

    corr_t = ak.concatenate([times, rc_data_ch.t0], axis=1)
    corr_a = ak.concatenate([amps, rc_data_ch.npe], axis=1)

    return corr_t, corr_a
