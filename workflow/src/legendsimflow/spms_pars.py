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
from pathlib import Path

import awkward as ak
import numpy as np
from lgdo import lh5

from .profile import make_profiler
from .utils import lookup_dataflow_config

log = logging.getLogger(__name__)


def lookup_evt_files(l200data: str, runid: str, evt_tier_name: str) -> list[str | Path]:
    """Lookup the paths to the `evt` files."""
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


def _next_rc_evt_file(evt_files: list[str | Path], rc_file_state: dict) -> str | Path:
    """Return the next `evt` file, cycling through files in input order before repeating."""
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
    rc_evt_files: list[str | Path],
) -> dict[str, dict[str, np.ndarray]]:
    """Build per-file trigger index lookup for RC extraction.

    Returns a dictionary keyed by file path string with entries:
    - ``forced_pulser``: indices for forced/pulser and non-muon events
    - ``geds``: indices for geds and non-muon events
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
    rc_evt_files: list[str | Path],
    rc_file_state: dict,
    chunk_size: int,
    evt_tier_name: str,
    hit_tier_name: str,
    sipm: str,
    sipm_uid: int,
    rc_index_lookup: dict[str, dict[str, np.ndarray]],
) -> ak.Array:
    """Assemble random-coincidence data for one chunk.

    Parameters
    ----------
    rc_evt_files
        Ordered list of evt files that can provide random-coincidence data.
    rc_file_state
        Mutable state for file cycling and carryover between chunks. Expected
        keys are created/updated internally (e.g. ``order``, ``idx``,
        ``counts``, ``carryover``).
    chunk_size
        Number of random-coincidence events requested for the current chunk.
    evt_tier_name
        Tier name of the evt file path.
    hit_tier_name
        Tier name used to derive the hit file path from the evt file path.
    sipm
        SiPM channel name. Use ``"all"`` for the summed SiPM stream.
    sipm_uid
        SiPM channel UID used when selecting a specific channel.
    rc_index_lookup
        Precomputed mapping from evt file to trigger-event indices,
        built with ``build_rc_evt_index_lookup``.

    Returns
    -------
    ak.Array
        Random-coincidence data for one chunk with fields ``npe`` and ``t0``.
    """
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
            evt_tier_name,
            evt_tier_name,
            hit_tier_name,
            sipm,
            sipm_uid,
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
    win_ranges: list[tuple[float, float]],
    time_domain_ns: tuple[float, float],
    min_sep_ns: float,
) -> tuple[ak.Array, ak.Array]:
    """Helper function to process SiPM data within specified window ranges.

    Parameters
    ----------
    time
        SiPM `t0` array from `evt` file, or equivalently, `trigger_pos` with `is_valid_hit` from `hit` file.
    energy
        SiPM `energy` array from `evt` file, or equivalently, `energy_in_pe` with `is_valid_hit` from `hit` file.
    win_ranges
        List of `(start, end)` tuples defining window ranges in nanoseconds.
    time_domain_ns
        Target time range `(start, end)` for output times in nanoseconds.
        E.g., `(-1000, 5000)` means output times will be in `[-1000, 5000]`.
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
            tmsk = (time >= wstart) & (time < wend)
            npe_tmp = energy[tmsk]
            t0_tmp = time[tmsk] - (wstart - time_domain_ns[0])

            npe_list.append(npe_tmp)
            t0_list.append(t0_tmp)

    if not npe_list:
        return ak.Array([]), ak.Array([])

    return ak.concatenate(npe_list), ak.concatenate(t0_list)


def get_rc_evt_mask(evt_file: str | Path) -> tuple[ak.Array, ak.Array]:
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
    evt_tier_name: str,
    hit_tier_name: str,
    sipm: str,
    sipm_uid: int,
    rc_index_lookup: dict[str, dict[str, np.ndarray]],
    time_domain_ns: tuple[float, float] = (-1_000, 5_000),
    min_sep_ns: float = 6_000,
    ext_trig_range_ns: list[tuple[float, float]] | None = None,
    ge_trig_range_ns: list[tuple[float, float]] | None = None,
) -> ak.Array:
    """Extract a library of forced trigger events.

    To be used in correcting the SiPM photoelectrons with random coincidences.

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
    evt_file
        Event tier data file.
    evt_tier_name
        Tier name of the evt file path.
    hit_tier_name
        Tier name used to derive the hit file path from the evt file path.
    sipm
        SiPM channel name to extract data for. If "all", flatten channel
        dimension.
    sipm_uid
        SiPM channel ID to extract data for.
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
        Window ranges for forced/pulser trigger events, as list of (start, end)
        tuples in nanoseconds. Default: ``[(1_000, 44_000), (55_000,
        100_000)]``.
    ge_trig_range_ns
        Window ranges for HPGe/LAr trigger events, as list of (start, end)
        tuples in nanoseconds.  Default: ``[(1_000, 44_000)]``.

    Returns
    -------
    Array with fields "npe", the number of pe per hit, and "t0", the
    corresponding time relative to the start of a window in the trace (bounded
    by ``time_domain_ns``).
    """
    perf_block, print_perf, _ = make_profiler()

    npe_list: list[ak.Array] = []
    t0_list: list[ak.Array] = []
    n_collected_events = 0

    # Set defaults if not provided
    if ext_trig_range_ns is None:
        ext_trig_range_ns = [
            (1_000, 44_000),
            (55_000, 100_000),
        ]  # forced/pulser events with full waveform windows except the central window around the trigger
    if ge_trig_range_ns is None:
        ge_trig_range_ns = [
            (1_000, 44_000)
        ]  # geds trigger events only in the window before the trigger

    evt_file_key = str(evt_file)
    if evt_file_key in rc_index_lookup:
        idx_fp = rc_index_lookup[evt_file_key]["forced_pulser"]
        idx_getrg = rc_index_lookup[evt_file_key]["geds"]
    else:
        # Fallback if key not found (shouldn't happen if build_rc_evt_index_lookup was complete)
        mask_fp, mask_getrg = get_rc_evt_mask(evt_file)
        idx_fp = ak.where(mask_fp)[0].to_numpy()
        idx_getrg = ak.where(mask_getrg)[0].to_numpy()

    n_forced_pulser = len(idx_fp)
    n_geds = len(idx_getrg)

    if sipm == "all":
        with perf_block("ftlib_read_evt_spms()"):
            evt = lh5.read(
                "evt",
                evt_file,
                field_mask=[
                    "spms/energy",
                    "spms/t0",
                ],
            ).view_as("ak")

            spms_fp = evt.spms[idx_fp]
            time_fp = ak.flatten(spms_fp.t0, axis=-1)
            energy_fp = ak.flatten(spms_fp.energy, axis=-1)

            spms_getrg = evt.spms[idx_getrg]
            time_getrg = ak.flatten(spms_getrg.t0, axis=-1)
            energy_getrg = ak.flatten(spms_getrg.energy, axis=-1)

    else:
        tcm_file = Path(str(evt_file).replace(evt_tier_name, "tcm"))
        hit_file = Path(str(evt_file).replace(evt_tier_name, hit_tier_name))

        with perf_block("ftlib_read_hit_table()"):
            data_ch = lh5.read_as(
                f"ch{sipm_uid}/hit/",
                hit_file,
                "ak",
                field_mask=["trigger_pos", "energy_in_pe", "is_valid_hit"],
            )

        with perf_block("ftlib_read_tcm_tables()"):
            idx_union = np.unique(np.concatenate((idx_fp, idx_getrg)))
            tcm_all = lh5.read_as(
                "hardware_tcm_1",
                tcm_file,
                "ak",
                idx=idx_union,
                field_mask=["table_key", "row_in_table"],
            )
            tcm_fp = tcm_all[np.searchsorted(idx_union, idx_fp)]
            tcm_getrg = tcm_all[np.searchsorted(idx_union, idx_getrg)]

        with perf_block("ftlib_extract_channel_rows()"):
            mask = tcm_fp.table_key == sipm_uid
            rows = ak.flatten(tcm_fp.row_in_table[mask]).to_numpy()

            if len(rows) > 0:
                data_ch_fp = data_ch[rows]
                time_fp = data_ch_fp.trigger_pos[data_ch_fp.is_valid_hit]
                energy_fp = data_ch_fp.energy_in_pe[data_ch_fp.is_valid_hit]
            else:
                time_fp = ak.Array([])
                energy_fp = ak.Array([])

            mask = tcm_getrg.table_key == sipm_uid
            rows = ak.flatten(tcm_getrg.row_in_table[mask]).to_numpy()

            if len(rows) > 0:
                data_ch_getrg = data_ch[rows]
                time_getrg = data_ch_getrg.trigger_pos[data_ch_getrg.is_valid_hit]
                energy_getrg = data_ch_getrg.energy_in_pe[data_ch_getrg.is_valid_hit]
            else:
                time_getrg = ak.Array([])
                energy_getrg = ak.Array([])

    if len(time_fp) > 0:
        with perf_block("ftlib_process_windows()"):
            npe_chunk, t0_chunk = _process_spms_windows(
                time_fp, energy_fp, ext_trig_range_ns, time_domain_ns, min_sep_ns
            )
            if len(npe_chunk) > 0:
                npe_list.append(npe_chunk)
                t0_list.append(t0_chunk)
                n_collected_events += len(npe_chunk)
    if len(time_getrg) > 0:
        with perf_block("ftlib_process_windows()"):
            npe_chunk, t0_chunk = _process_spms_windows(
                time_getrg, energy_getrg, ge_trig_range_ns, time_domain_ns, min_sep_ns
            )
            if len(npe_chunk) > 0:
                npe_list.append(npe_chunk)
                t0_list.append(t0_chunk)
                n_collected_events += len(npe_chunk)

    log.debug(
        "forced-trigger library file %s: forced_or_pulser_events=%d "
        "geds_events=%d cumulative_events=%d",
        Path(evt_file).name,
        n_forced_pulser,
        n_geds,
        n_collected_events,
    )

    with perf_block("ftlib_concatenate()"):
        if npe_list:
            npe = ak.concatenate(npe_list)
            t0 = ak.concatenate(t0_list)
        else:
            npe = ak.Array([])
            t0 = ak.Array([])

    print_perf()

    return ak.Array({"npe": npe, "t0": t0})
