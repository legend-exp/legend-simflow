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
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import awkward as ak
import numpy as np
from lgdo import lh5

from .profile import make_profiler
from .utils import lookup_dataflow_config

log = logging.getLogger(__name__)


def lookup_evt_files(
    l200data: str | Path, runid: str, evt_tier_name: str
) -> list[Path]:
    """Look up the `evt` tier file paths for a given run.

    Parameters
    ----------
    l200data
        Root path to the LEGEND-200 data directory.
    runid
        Run identifier string (e.g. ``"l200-p16-r008-phy"``).
    evt_tier_name
        Name of the evt tier (e.g. ``"evt"``).

    Returns
    -------
    list[Path]
        Matching evt-tier file paths for the given run.
    """
    _, period, run, data_type = re.split(r"\W+", runid)

    if isinstance(l200data, str):
        l200data = Path(l200data)

    dataflow_config = lookup_dataflow_config(l200data)

    # get the paths to evt tier files
    df_cfg = (
        dataflow_config["setups"]["l200"]["paths"]
        if ("setups" in dataflow_config)
        else dataflow_config["paths"]
    )

    evt_path = Path(df_cfg[f"tier_{evt_tier_name}"]).resolve()
    return list((evt_path / data_type / period / run).glob("*"))


def _next_rc_evt_file(
    evt_files: Sequence[str | Path], rc_file_state: dict[str, Any]
) -> str | Path:
    """Return the next evt file, cycling through the list before repeating.

    Parameters
    ----------
    evt_files
        Ordered sequence of evt file paths to cycle through.
    rc_file_state
        Mutable state dict shared across calls.  On the first call it is
        populated with keys ``order`` (the file list), ``idx`` (current
        position, int), and ``completed_cycle`` (bool, set to ``True`` once
        the list has been exhausted once).  Subsequent calls increment ``idx``
        and wrap it when all files have been visited.

    Returns
    -------
    str | Path
        Path to the next evt file to process.
    """
    if "order" not in rc_file_state:
        order = list(evt_files)
        rc_file_state["order"] = order
        rc_file_state["idx"] = 0
        rc_file_state["completed_cycle"] = False

    order = rc_file_state["order"]
    idx = rc_file_state["idx"]

    # if all files were used once, start again from the first file
    if idx >= len(order):
        idx = 0
        if not rc_file_state["completed_cycle"]:
            log.warning(
                "restarting evt file iteration for RC correction; "
                "cycling through %d files again",
                len(order),
            )
            rc_file_state["completed_cycle"] = True

    evt_file = order[idx]
    rc_file_state["idx"] = idx + 1

    return evt_file


def build_rc_evt_index_lookup(
    rc_evt_files: Sequence[str | Path],
) -> dict[str, dict[str, np.ndarray]]:
    """Build per-file trigger index lookup for RC extraction.

    Parameters
    ----------
    rc_evt_files
        Evt-tier files to index.

    Returns
    -------
    dict
        Dictionary keyed by file path (as string) with entries:

        - ``forced_pulser``: row indices of forced/pulser, non-muon events
        - ``geds``: row indices of HPGe-triggered, non-muon events
    """
    lookup: dict[str, dict[str, np.ndarray]] = {}
    for evt_file in rc_evt_files:
        mask_fp, mask_getrg = get_rc_evt_mask(evt_file)
        lookup[str(evt_file)] = {
            "forced_pulser": ak.where(mask_fp)[0].to_numpy(),
            "geds": ak.where(mask_getrg)[0].to_numpy(),
        }
    return lookup


def get_chunk_rc_data(
    rc_evt_files: Sequence[str | Path],
    rc_file_state: dict[str, Any],
    chunk_size: int,
    rc_index_lookup: dict[str, dict[str, np.ndarray]],
) -> ak.Array:
    """Assemble random-coincidence data for one chunk.

    Parameters
    ----------
    rc_evt_files
        Ordered sequence of evt files that can provide random-coincidence data.
        Must not be empty.
    rc_file_state
        Mutable state for file cycling and carryover between chunks. Expected
        keys are created/updated internally (e.g. ``order``, ``idx``,
        ``counts``, ``carryover``).
    chunk_size
        Number of random-coincidence events requested for the current chunk.
        Must be positive.
    rc_index_lookup
        Precomputed mapping from evt file to trigger-event indices,
        built with :func:`build_rc_evt_index_lookup`.

    Returns
    -------
    ak.Array
        Random-coincidence data for one chunk with fields ``rawid``
        ``(chunk_size, n_channels)``, ``npe`` ``(chunk_size, n_channels, n_pe)``
        and ``t0`` (same shape as ``npe``).
    """
    if not rc_evt_files:
        msg = "rc_evt_files must not be empty"
        raise ValueError(msg)
    if chunk_size <= 0:
        msg = f"chunk_size must be positive, got {chunk_size}"
        raise ValueError(msg)

    rc_parts: list[ak.Array] = []
    total_rc_events = 0
    empty_parts_streak = 0
    max_empty_parts = max(2 * len(rc_evt_files), 1)

    carryover = rc_file_state.get("carryover")
    if carryover is not None and len(carryover) > 0:
        carryover_len = len(carryover)
        n_from_carryover = min(chunk_size, len(carryover))
        rc_parts.append(carryover[:n_from_carryover])
        total_rc_events += n_from_carryover
        if n_from_carryover < len(carryover):
            rc_file_state["carryover"] = carryover[n_from_carryover:]
        else:
            rc_file_state["carryover"] = None
        remaining_carryover = rc_file_state.get("carryover")
        remaining_carryover_len = (
            len(remaining_carryover) if remaining_carryover is not None else 0
        )
        log.debug(
            "used %d residual RC events from carryover; %d remained queued "
            "(had %d, chunk needs %d total)",
            n_from_carryover,
            remaining_carryover_len,
            carryover_len,
            chunk_size,
        )

    while total_rc_events < chunk_size and empty_parts_streak < max_empty_parts:
        rc_evt_file = _next_rc_evt_file(rc_evt_files, rc_file_state)
        n_missing = chunk_size - total_rc_events
        part = get_rc_library(
            rc_evt_file,
            rc_index_lookup,
        )
        if len(part) == 0:
            empty_parts_streak += 1
            log.warning(
                "forced-trigger library from %s is empty for this chunk; "
                "trying next file (%d/%d consecutive empties)",
                Path(rc_evt_file).name,
                empty_parts_streak,
                max_empty_parts,
            )
            continue
        empty_parts_streak = 0

        n_take = min(n_missing, len(part))
        rc_parts.append(part[:n_take])
        total_rc_events += n_take

        n_left = len(part) - n_take
        if n_left > 0:
            rc_file_state["carryover"] = part[n_take:]
            log.debug(
                "stored %d residual RC events from %s after taking %d/%d "
                "for this chunk",
                n_left,
                Path(rc_evt_file).name,
                n_take,
                len(part),
            )

    if total_rc_events == 0:
        msg = "no random coincidences available from any evt file"
        raise RuntimeError(msg)
    if total_rc_events < chunk_size:
        msg = (
            "insufficient random coincidences to fill chunk: "
            f"needed {chunk_size}, obtained {total_rc_events} "
            f"after {empty_parts_streak} consecutive empty libraries"
        )
        raise RuntimeError(msg)

    return ak.concatenate(rc_parts) if len(rc_parts) > 1 else rc_parts[0]


def _process_spms_windows(
    time: ak.Array,
    energy: ak.Array,
    win_ranges: Sequence[tuple[float, float]],
    time_domain_ns: tuple[float, float],
    min_sep_ns: float,
) -> tuple[ak.Array, ak.Array]:
    """Process SiPM data within specified window ranges.

    Each ``(start, end)`` range in ``win_ranges`` is tiled with non-overlapping
    windows of length ``time_domain_ns[1] - time_domain_ns[0]``, separated by
    ``min_sep_ns``.  PE hits falling inside each window are selected and their
    times are shifted so that the window start maps to ``time_domain_ns[0]``.

    The function works on arrays of any rank.  For N-D input (e.g. shape
    ``(n_events, n_channels, n_pe)``), each extracted window produces one
    output block of the same shape along all but the innermost axis, with only
    the PE dimension filtered.  Blocks from all windows are then concatenated
    along axis=0, so M source events processed through W windows yield
    ``W * M`` output entries.

    Parameters
    ----------
    time
        PE hit times.  Any shape; the innermost axis is the PE axis.
    energy
        PE energies, same shape as ``time``.
    win_ranges
        List of ``(start, end)`` tuples defining the time ranges to tile, in
        nanoseconds.
    time_domain_ns
        Target time range ``(start, end)`` for output times in nanoseconds.
        The window length is ``end - start``.  E.g. ``(-1000, 5000)`` selects
        6000 ns windows and maps their start to ``-1000 ns``.
    min_sep_ns
        Minimum gap between consecutive windows in nanoseconds.

    Returns
    -------
    npe
        PE energies extracted from all windows, concatenated along axis=0.
    t0
        PE times relative to each window's start (bounded by
        ``time_domain_ns``), same shape as ``npe``.
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
            tmsk = (time >= wstart) & (time < wend)
            npe_tmp = energy[tmsk]
            t0_tmp = time[tmsk] - (wstart - time_domain_ns[0])

            npe_list.append(npe_tmp)
            t0_list.append(t0_tmp)

    if not npe_list:
        return ak.Array([]), ak.Array([])

    return ak.concatenate(npe_list), ak.concatenate(t0_list)


def get_rc_evt_mask(evt_file: str | Path) -> tuple[ak.Array, ak.Array]:
    """Compute boolean event masks for random-coincidence extraction.

    Parameters
    ----------
    evt_file
        Path to the evt-tier LH5 file.

    Returns
    -------
    mask_forced_pulser
        Boolean mask selecting forced-trigger and pulser events, excluding
        muon coincidences.
    mask_geds
        Boolean mask selecting HPGe-triggered events, excluding muon
        coincidences.
    """
    evt = lh5.read(
        "evt",
        evt_file,
        field_mask=[
            "trigger/is_forced",
            "coincident/geds",
            "coincident/muon_offline",
            "coincident/puls",
        ],
    ).view_as("ak")

    is_forced = evt.trigger.is_forced
    is_geds_trig = evt.coincident.geds
    is_muon = evt.coincident.muon_offline
    is_pulser = evt.coincident.puls  # codespell:ignore puls

    mask_forced_pulser = (is_forced | is_pulser) & ~is_muon
    mask_geds = is_geds_trig & ~is_muon

    return mask_forced_pulser, mask_geds


def get_rc_library(
    evt_file: str | Path,
    rc_index_lookup: dict[str, dict[str, np.ndarray]],
    time_domain_ns: tuple[float, float] = (-1_000, 5_000),
    min_sep_ns: float = 6_000,
    ext_trig_range_ns: Sequence[tuple[float, float]] = (
        (1_000, 44_000),
        (55_000, 100_000),
    ),
    ge_trig_range_ns: Sequence[tuple[float, float]] = ((1_000, 44_000),),
) -> ak.Array:
    """Extract a library of random-coincidence (RC) events from an evt file.

    To be used in correcting the SiPM photoelectrons with random coincidences.

    For each qualifying trigger event, the SiPM waveform is divided into
    multiple non-overlapping time windows (see ``_process_spms_windows``).
    Each window yields one independent RC event, so the total number of entries
    in the returned library is ``n_source_events x n_windows``.  The
    per-channel structure is preserved: ``npe`` and ``t0`` have shape
    ``(n_rc_events, n_channels, n_pe)`` and ``rawid`` has shape
    ``(n_rc_events, n_channels)``, matching the ``spms/*`` layout of the evt
    tier.

    Two trigger categories are processed with different window ranges to avoid
    contaminating RC events with physics signal:

    - Forced/pulser triggers: full waveform outside the central trigger window,
      ``((1_000, 44_000), (55_000, 100_000))`` ns by default.
    - HPGe/LAr triggers: first half only (before the trigger), ``((1_000,
      44_000),)`` ns by default.

    Both categories are filtered to exclude muon coincidences.

    Parameters
    ----------
    evt_file
        Event tier data file.
    rc_index_lookup
        Precomputed mapping from evt file to trigger-event indices,
        built with ``build_rc_evt_index_lookup``.
    time_domain_ns
        Target time range (start, end) for output times in nanoseconds.  E.g.,
        ``(-1000, 5000)`` means output times will be in ``[-1000, 5000]``.
        Default: ``(-1_000, 5_000)``.
    min_sep_ns
        Minimal separation time between two windows in a trace, in nanoseconds.
        Default 6000.
    ext_trig_range_ns
        Window ranges for forced/pulser trigger events, as a sequence of
        ``(start, end)`` pairs in nanoseconds. Default:
        ``((1_000, 44_000), (55_000, 100_000))``.
    ge_trig_range_ns
        Window ranges for HPGe/LAr trigger events, as a sequence of
        ``(start, end)`` pairs in nanoseconds. Default: ``((1_000, 44_000),)``.

    Returns
    -------
    ak.Array
        Record array with fields ``rawid`` (channel UIDs, shape
        ``(n_rc_events, n_channels)``), ``npe`` (PE energies, shape
        ``(n_rc_events, n_channels, n_pe)``), and ``t0`` (times relative to
        each window start, same shape as ``npe``).  Channel ordering within
        each event matches the source ``spms/rawid`` ordering in ``evt_file``.
    """
    perf_block, _, _ = make_profiler()

    evt_file_key = str(evt_file)
    idx_fp = rc_index_lookup[evt_file_key]["forced_pulser"]
    idx_getrg = rc_index_lookup[evt_file_key]["geds"]

    n_forced_pulser = len(idx_fp)
    n_geds = len(idx_getrg)

    if n_forced_pulser == 0 and n_geds == 0:
        log.debug("no forced/pulser or geds events found in %s", Path(evt_file).name)
        return ak.Array(
            {"rawid": ak.Array([]), "npe": ak.Array([]), "t0": ak.Array([])}
        )

    with perf_block("ftlib_read_evt_spms()"):
        evt = lh5.read(
            "evt",
            evt_file,
            field_mask=["spms/energy", "spms/t0", "spms/rawid"],
        ).view_as("ak")

    results: list[ak.Array] = []
    n_collected_events = 0

    for idx, win_ranges in [
        (idx_fp, ext_trig_range_ns),
        (idx_getrg, ge_trig_range_ns),
    ]:
        if len(idx) == 0 or not win_ranges:
            continue
        spms = evt.spms[idx]
        with perf_block("ftlib_process_windows()"):
            npe, t0 = _process_spms_windows(
                spms.t0, spms.energy, win_ranges, time_domain_ns, min_sep_ns
            )
        if len(npe) > 0:
            # each window produced one RC event per source event; repeat rawid
            # accordingly so it aligns with the (n_windows*n_src, n_ch, n_pe) output
            n_windows = len(npe) // len(idx)
            rawid = ak.concatenate([spms.rawid] * n_windows)
            results.append(ak.Array({"rawid": rawid, "npe": npe, "t0": t0}))
            n_collected_events += len(npe)

    log.debug(
        "forced-trigger library file %s: forced_or_pulser_events=%d "
        "geds_events=%d cumulative_events=%d",
        Path(evt_file).name,
        n_forced_pulser,
        n_geds,
        n_collected_events,
    )

    if not results:
        return ak.Array(
            {"rawid": ak.Array([]), "npe": ak.Array([]), "t0": ak.Array([])}
        )

    return ak.concatenate(results) if len(results) > 1 else results[0]
