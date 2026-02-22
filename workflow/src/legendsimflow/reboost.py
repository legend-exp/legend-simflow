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
import random
import re
import time
from collections.abc import Iterable, Mapping
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


def _process_spms_windows(
    spms: ak.Array,
    win_ranges: list[tuple[float, float]],
    time_domain_ns: tuple[float, float],
    min_sep_ns: float,
) -> tuple[ak.Array, ak.Array]:
    """Helper function to process SiPM data within specified window ranges.

    Parameters
    ----------
    spms
        SiPM data array with fields 'energy' and 't0'.
    win_ranges
        List of (start, end) tuples defining window ranges in nanoseconds.
    time_domain_ns
        Target time range (start, end) for output times in nanoseconds.
        E.g., (-1000, 5000) means output times will be in [-1000, 5000].
    min_sep_ns
        Minimal separation between windows in nanoseconds.
    Returns
    -------
    npe
        Photoelectron counts extracted from the requested windows.
    t0
        Times relative to time_domain_ns extracted from the requested windows.
    """
    # Validate inputs to avoid infinite loops and invalid window definitions
    if time_domain_ns[1] <= time_domain_ns[0]:
        msg = (
            f"time_domain_ns must have time_domain_ns[1] > time_domain_ns[0], "
            f"got {time_domain_ns}"
        )
        raise ValueError(msg)
    if min_sep_ns < 0:
        msg = f"min_sep_ns must be non-negative, got {min_sep_ns}"
        raise ValueError(msg)
    win_len_ns = time_domain_ns[1] - time_domain_ns[0]

    npe_list = []
    t0_list = []

    for win_range in win_ranges:
        # Ensure each window range is a valid (start, end) pair
        if len(win_range) != 2:
            msg = f"Each win_range must be a (start, end) tuple, got {win_range!r}"
            raise ValueError(msg)
        start, end = win_range
        if end <= start:
            msg = (
                f"Each win_range must satisfy start < end, got start={start}, end={end}"
            )
            raise ValueError(msg)
        starts = []
        current_start = start
        while current_start + win_len_ns <= end:
            starts.append(current_start)
            current_start += win_len_ns + min_sep_ns
        ends = [s + win_len_ns for s in starts]

        for wstart, wend in zip(starts, ends, strict=True):
            tmsk = (spms.t0 >= wstart) & (spms.t0 < wend)
            npe_tmp = spms.energy[tmsk]
            t0_tmp = spms.t0[tmsk] - (wstart - time_domain_ns[0])

            npe_list.append(npe_tmp)
            t0_list.append(t0_tmp)

    if not npe_list:
        return ak.Array([]), ak.Array([])

    return ak.concatenate(npe_list), ak.concatenate(t0_list)


def get_forced_trigger_library(
    evt_files: Iterable[str],
    min_num_evts: int,
    time_domain_ns: tuple[float, float] = (-1_000, 5_000),
    min_sep_ns: float = 6_000,
    ext_trig_range_ns: list[tuple[float, float]] | None = None,
    ge_trig_range_ns: list[tuple[float, float]] | None = None,
) -> ak.Array:
    """Extract a library of forced trigger events to be used
    in correcting the SiPM pe and times with random coincidences.

    This reformats the data to make use of several windows within a waveform to
    build forced trigger events and then stores the number of pe and times for
    each SiPM channel from that window (to be used for corrections).

    This function processes two types of triggers with different window ranges:
    - Forced/pulser triggers: uses full waveform ``[(1_000, 44_000), (55_000,
      100_000)]`` ns
    - HPGe/LAr triggers: uses first half only ``(1_000, 44_000)`` ns

    Both are always filtered to exclude muon coincidences.

    Parameters
    ----------
    evt_files
        List of event tier data files.
    min_num_evts:
        Number of events required for forced trigger correction.
        Function will return at least ``min_num_evts``.
    time_domain_ns
        Target time range (start, end) for output times in nanoseconds.  E.g.,
        ``(-1000, 5000)`` means output times will be in ``[-1000, 5000]``.
        Default: ``(-1_000, 5_000)``.
    min_sep_ns
        Minimal separation time between two windows in a trace, in nanoseconds.
    ext_trig_range_ns
        Window ranges for forced/pulser trigger events, as list of (start, end)
        tuples in nanoseconds. Default: ``[(1_000, 44_000), (55_000,
        100_000)]``.
    ge_trig_range_ns
        Window ranges for HPGe/LAr trigger events, as list of (start, end)
        tuples in nanoseconds.  Default: ``[(1_000, 44_000)]``.

    Returns
    -------
    Array with fields "npe", the number of pe per SiPM and per hit, "t0", the
    time relative to the start of a window in the trace, per SiPM and per hit
    (makes sure that t0 are between bounds specified in time_domain_ns), and
    "rawid" the SiPM channel numbers.

    Example
    -------
    rc_data = get_forced_trigger_library(evt_files)
    """

    npe_chunks: list[ak.Array] = []
    t0_chunks: list[ak.Array] = []
    n_collected_events = 0
    rawids = None

    # Set defaults if not provided
    if ext_trig_range_ns is None:
        ext_trig_range_ns = [(1_000, 44_000), (55_000, 100_000)]
    if ge_trig_range_ns is None:
        ge_trig_range_ns = [(1_000, 44_000)]

    # Shuffle evt_files and process until we have enough events
    evt_files = list(evt_files)
    random.shuffle(evt_files)

    t_start_total = time.perf_counter()
    files_processed = 0
    total_read_wall_s = 0.0
    total_forced_pulser_wall_s = 0.0
    total_geds_wall_s = 0.0

    for file in evt_files:
        if n_collected_events >= min_num_evts:
            break

        files_processed += 1

        # Load all necessary data once
        t_read_start = time.perf_counter()
        evt = lh5.read(
            "evt/",
            file,
            field_mask=[
                "trigger/is_forced",
                "coincident/geds",
                "coincident/muon_offline",
                "coincident/puls",  # codespell:ignore puls
                "spms/energy",
                "spms/t0",
                "spms/rawid",
            ],
        ).view_as("ak")
        file_read_wall_s = time.perf_counter() - t_read_start
        total_read_wall_s += file_read_wall_s

        is_forced = evt.trigger.is_forced
        is_geds_trig = evt.coincident.geds
        is_muon = evt.coincident.muon_offline
        is_pulser = evt.coincident.puls  # codespell:ignore puls

        rawids_tmp = evt.spms.rawid[0]

        if rawids is not None and not ak.all(rawids == rawids_tmp):
            msg = "rawid should be the same in all cases"
            raise ValueError(msg)

        # Process forced/pulser events with full waveform windows
        mask_forced_pulser = (is_forced | is_pulser) & ~is_muon
        n_forced_pulser = int(ak.sum(mask_forced_pulser))
        spms_fp = evt.spms[mask_forced_pulser]
        file_forced_pulser_wall_s = 0.0
        if len(spms_fp) > 0:
            t_forced_start = time.perf_counter()
            npe_chunk, t0_chunk = _process_spms_windows(
                spms_fp, ext_trig_range_ns, time_domain_ns, min_sep_ns
            )
            if len(npe_chunk) > 0:
                npe_chunks.append(npe_chunk)
                t0_chunks.append(t0_chunk)
                n_collected_events += len(npe_chunk)
            file_forced_pulser_wall_s = time.perf_counter() - t_forced_start
            total_forced_pulser_wall_s += file_forced_pulser_wall_s
        # Process geds trigger events with limited window
        mask_geds = is_geds_trig & ~is_muon
        n_geds = int(ak.sum(mask_geds))
        spms_ge_trig = evt.spms[mask_geds]
        file_geds_wall_s = 0.0
        if len(spms_ge_trig) > 0:
            t_geds_start = time.perf_counter()
            npe_chunk, t0_chunk = _process_spms_windows(
                spms_ge_trig, ge_trig_range_ns, time_domain_ns, min_sep_ns
            )
            if len(npe_chunk) > 0:
                npe_chunks.append(npe_chunk)
                t0_chunks.append(t0_chunk)
                n_collected_events += len(npe_chunk)
            file_geds_wall_s = time.perf_counter() - t_geds_start
            total_geds_wall_s += file_geds_wall_s

        log.debug(
            "forced-trigger library file %s: read_wall_s=%.4g forced_or_pulser_events=%d "
            "forced_or_pulser_wall_s=%.4g geds_events=%d geds_wall_s=%.4g cumulative_events=%d",
            Path(file).name,
            file_read_wall_s,
            n_forced_pulser,
            file_forced_pulser_wall_s,
            n_geds,
            file_geds_wall_s,
            n_collected_events,
        )

        rawids = rawids_tmp

    total_wall_s = time.perf_counter() - t_start_total
    log.debug(
        "forced-trigger library summary: files_processed=%d total_wall_s=%.4g read_wall_s=%.4g "
        "forced_or_pulser_wall_s=%.4g geds_wall_s=%.4g returned_events=%d requested_min_events=%d",
        files_processed,
        total_wall_s,
        total_read_wall_s,
        total_forced_pulser_wall_s,
        total_geds_wall_s,
        n_collected_events,
        min_num_evts,
    )

    if npe_chunks:
        npe = ak.concatenate(npe_chunks)
        t0 = ak.concatenate(t0_chunks)
    else:
        npe = ak.Array([])
        t0 = ak.Array([])

    # Handle case where no events passed the filters
    if len(npe) == 0 or rawids is None:
        log.warning(
            "No events passed the filters in get_forced_trigger_library, returning empty arrays"
        )
        rawid = np.empty((0, len(rawids) if rawids is not None else 0), dtype=np.int32)
    else:
        rawid = np.vstack([rawids] * len(npe))

    return ak.Array(
        {
            "npe": npe,
            "t0": t0,
            "rawid": rawid,
        }
    )
