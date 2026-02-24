from __future__ import annotations

import logging
import random
import re
from collections.abc import Iterable
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
        SiPM data array with fields `energy` and `t0`.
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
            tmsk = (spms.t0 >= wstart) & (spms.t0 < wend)
            npe_tmp = spms.energy[tmsk]
            t0_tmp = spms.t0[tmsk] - (wstart - time_domain_ns[0])

            npe_list.append(npe_tmp)
            t0_list.append(t0_tmp)

    if not npe_list:
        return ak.Array([]), ak.Array([])

    return ak.concatenate(npe_list), ak.concatenate(t0_list)


def get_rc_library(
    evt_files: Iterable[str],
    min_num_evts: int,
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
    evt_files
        List of event tier data files.
    min_num_evts
        Minimum number of events requested for forced trigger correction. The
        function attempts to collect at least ``min_num_evts`` events but may
        return fewer if insufficient data are available in the input files.
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
    Array with fields "npe", the number of pe per SiPM and per hit, "t0", the
    time relative to the start of a window in the trace, per SiPM and per hit
    (makes sure that t0 are between bounds specified in `time_domain_ns`), and
    "rawid" the SiPM channel numbers.
    """
    perf_block, print_perf = make_profiler()

    npe_chunks: list[ak.Array] = []
    t0_chunks: list[ak.Array] = []
    n_collected_events = 0
    rawids = None

    # Set defaults if not provided
    if ext_trig_range_ns is None:
        ext_trig_range_ns = [(1_000, 44_000), (55_000, 100_000)]
    if ge_trig_range_ns is None:
        ge_trig_range_ns = [(1_000, 44_000)]

    # shuffle evt_files in case rc change during the run
    evt_files = list(evt_files)
    random.shuffle(evt_files)

    files_processed = 0

    for file in evt_files:
        if n_collected_events >= min_num_evts:
            break

        files_processed += 1

        # Load all necessary data once
        with perf_block("ftlib_read_data"):
            evt = lh5.read(
                "evt",
                file,
                field_mask=[
                    "trigger/is_forced",
                    "coincident/geds",
                    "coincident/muon_offline",
                    "coincident/puls",
                    "spms/energy",
                    "spms/t0",
                    "spms/rawid",
                ],
            ).view_as("ak")

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
        if len(spms_fp) > 0:
            with perf_block("ftlib_process_windows()"):
                npe_chunk, t0_chunk = _process_spms_windows(
                    spms_fp, ext_trig_range_ns, time_domain_ns, min_sep_ns
                )
                if len(npe_chunk) > 0:
                    npe_chunks.append(npe_chunk)
                    t0_chunks.append(t0_chunk)
                    n_collected_events += len(npe_chunk)

        # Process geds trigger events with limited window
        mask_geds = is_geds_trig & ~is_muon
        n_geds = int(ak.sum(mask_geds))
        spms_ge_trig = evt.spms[mask_geds]
        if len(spms_ge_trig) > 0:
            with perf_block("ftlib_process_windows()"):
                npe_chunk, t0_chunk = _process_spms_windows(
                    spms_ge_trig, ge_trig_range_ns, time_domain_ns, min_sep_ns
                )
                if len(npe_chunk) > 0:
                    npe_chunks.append(npe_chunk)
                    t0_chunks.append(t0_chunk)
                    n_collected_events += len(npe_chunk)

        log.debug(
            "forced-trigger library file %s: forced_or_pulser_events=%d "
            "geds_events=%d cumulative_events=%d",
            Path(file).name,
            n_forced_pulser,
            n_geds,
            n_collected_events,
        )

        rawids = rawids_tmp

    log.debug(
        "forced-trigger library summary: files_processed=%d "
        "returned_events=%d requested_min_events=%d",
        files_processed,
        n_collected_events,
        min_num_evts,
    )

    with perf_block("ftlib_concatenate()"):
        if npe_chunks:
            npe = ak.concatenate(npe_chunks)
            t0 = ak.concatenate(t0_chunks)
        else:
            npe = ak.Array([])
            t0 = ak.Array([])

        # Handle case where no events passed the filters
        if len(npe) == 0 or rawids is None:
            log.warning(
                "No events passed the filters in get_random_coincidences_library, "
                "returning empty arrays"
            )
            rawid = np.empty(
                (0, len(rawids) if rawids is not None else 0), dtype=np.int32
            )
        else:
            rawid = np.vstack([rawids] * len(npe))

    print_perf()

    return ak.Array(
        {
            "npe": npe,
            "t0": t0,
            "rawid": rawid,
        }
    )


def get_rc_library_chunk(
    rc_library: ak.Array,
    chunk_len: int,
    rc_offset: dict,
) -> ak.Array:
    """Select a chunk-length slice from a pre-sampled random coincidence library.

    Always returns exactly ``chunk_len`` entries using wrap-around indexing when
    necessary, and advances ``rc_offset["idx"]`` accordingly.
    """
    lib_len = len(rc_library)

    if chunk_len > lib_len:
        msg = (
            "forced trigger library smaller than chunk; "
            "reusing events with wrap-around."
        )
        log.warning(msg)

    start_idx = rc_offset.get("idx", 0)
    idx_array = (np.arange(chunk_len) + start_idx) % lib_len
    rc_offset["idx"] = int((start_idx + chunk_len) % lib_len)

    return rc_library[idx_array]


def get_sipm_rc_data(
    rc_library: ak.Array,
    sipm: str,
    sipm_uid: int,
) -> tuple[ak.Array, ak.Array]:
    """Extract data for a specific SiPM channel from the library.

    Filters the forced trigger library to return only the photoelectron counts
    and times for a single SiPM channel across all events.

    Parameters
    ----------
    ft_library
        Library of forced trigger events containing npe, t0, and rawid fields.
    sipm
        SiPM channel name to extract data for. If "all", flatten channel
        dimension.
    sipm_uid
        SiPM channel ID to extract data for.

    Returns
    -------
    npe
        Photoelectron counts for the SiPM channel across all events.
    t0
        Photoelectron times for the SiPM channel across all events.

    Raises
    ------
    ValueError
        If the SiPM UID is not found in the library.
    """
    if sipm == "all":
        npe = ak.flatten(rc_library.npe, axis=-1)
        t0 = ak.flatten(rc_library.t0, axis=-1)

    else:
        # Find the channel index for this SiPM UID
        # rawid[0] gives the channel IDs (should be the same for all events)
        channel_indices = ak.where(rc_library.rawid[0] == sipm_uid)[0]

        if len(channel_indices) == 0:
            msg = f"SiPM UID {sipm_uid} not found in forced trigger library"
            raise ValueError(msg)

        ch_idx = int(channel_indices[0])

        # Select data for this channel from all events
        npe = rc_library.npe[:, ch_idx]
        t0 = rc_library.t0[:, ch_idx]

    return npe, t0
