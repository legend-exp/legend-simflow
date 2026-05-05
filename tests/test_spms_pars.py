"""Tests for spms_pars module."""

from __future__ import annotations

import awkward as ak
import numpy as np
import pytest

from legendsimflow import spms_pars

EVT_FILE = "lh5/l200-p16-r008-ssc-20251006T205904Z-tier_evt.lh5"


def test_process_spms_windows_two_windows():
    """Verify window extraction and time shifting for a simple input."""
    time = ak.Array([2000, 3000, 7000, 8000, 9000])
    energy = ak.Array([1.0, 2.0, 5.0, 3.0, 4.0])

    npe, t0 = spms_pars._process_spms_windows(
        time,
        energy,
        [(1000, 14000)],
        (-1000, 5000),
        1000,
    )

    assert len(npe) == 4
    assert float(ak.sum(npe)) == 10.0
    flat_t0 = t0
    assert ak.all(flat_t0 >= -1000)
    assert ak.all(flat_t0 <= 5000)


def test_process_spms_windows_invalid_time_domain_raises():
    """time_domain_ns must be increasing."""
    with pytest.raises(ValueError, match="time_domain_ns"):
        spms_pars._process_spms_windows(
            ak.Array([1.0]),
            ak.Array([1.0]),
            [(0.0, 10.0)],
            (5.0, 5.0),
            0.0,
        )


def test_process_spms_windows_negative_min_sep_raises():
    """min_sep_ns must be non-negative."""
    with pytest.raises(ValueError, match="min_sep_ns"):
        spms_pars._process_spms_windows(
            ak.Array([1.0]),
            ak.Array([1.0]),
            [(0.0, 10.0)],
            (0.0, 5.0),
            -1.0,
        )


def test_process_spms_windows_invalid_range_raises():
    """Each window range must satisfy start < end."""
    with pytest.raises(ValueError, match="start < end"):
        spms_pars._process_spms_windows(
            ak.Array([1.0]),
            ak.Array([1.0]),
            [(10.0, 10.0)],
            (0.0, 5.0),
            0.0,
        )


def test_next_rc_evt_file_cycles_order():
    """_next_rc_evt_file should iterate in order and wrap around."""
    files = ["f0", "f1"]
    state: dict = {}

    assert spms_pars._next_rc_evt_file(files, state) == "f0"
    assert spms_pars._next_rc_evt_file(files, state) == "f1"
    assert spms_pars._next_rc_evt_file(files, state) == "f0"


def test_next_rc_evt_file_sets_cycle_flag_on_wrap():
    """completed_cycle flag is set when iteration wraps the first time."""
    files = ["f0"]
    state: dict = {}

    _ = spms_pars._next_rc_evt_file(files, state)
    assert state["completed_cycle"] is False
    _ = spms_pars._next_rc_evt_file(files, state)
    assert state["completed_cycle"] is True


# --- tests requiring a real evt file from legend-testdata ---


def test_get_rc_evt_mask(legend_testdata):
    evt_file = legend_testdata[EVT_FILE]
    mask_fp, mask_getrg = spms_pars.get_rc_evt_mask(evt_file)

    assert len(mask_fp) == 1000
    assert len(mask_getrg) == 1000
    # forced/pulser and geds triggers should be mutually exclusive categories
    # (forced events are never geds-triggered in this data)
    assert int(ak.sum(mask_fp)) > 0
    assert int(ak.sum(mask_getrg)) > 0


def test_build_rc_evt_index_lookup(legend_testdata):
    evt_file = legend_testdata[EVT_FILE]
    lookup = spms_pars.build_rc_evt_index_lookup([evt_file])

    assert str(evt_file) in lookup
    entry = lookup[str(evt_file)]
    assert "forced_pulser" in entry
    assert "geds" in entry
    assert isinstance(entry["forced_pulser"], np.ndarray)
    assert isinstance(entry["geds"], np.ndarray)
    assert len(entry["forced_pulser"]) > 0
    assert len(entry["geds"]) > 0


def test_get_rc_library(legend_testdata):
    evt_file = legend_testdata[EVT_FILE]
    lookup = spms_pars.build_rc_evt_index_lookup([evt_file])
    lib = spms_pars.get_rc_library(evt_file, lookup)

    assert "rawid" in lib.fields
    assert "npe" in lib.fields
    assert "t0" in lib.fields
    assert len(lib) > 0
    # rawid, npe, t0 share the same outer two dimensions (events x channels)
    assert len(lib.rawid) == len(lib.npe) == len(lib.t0)
    assert ak.all(lib.t0 >= -1000)
    assert ak.all(lib.t0 <= 5000)


def test_get_rc_library_window_multiplication(legend_testdata):
    """Each source event produces one RC event per extracted time window."""
    evt_file = legend_testdata[EVT_FILE]
    lookup = spms_pars.build_rc_evt_index_lookup([evt_file])
    n_fp = len(lookup[str(evt_file)]["forced_pulser"])

    # one window: exactly n_fp RC events
    result_1w = spms_pars.get_rc_library(
        evt_file,
        lookup,
        ext_trig_range_ns=[(1_000, 7_000)],
        ge_trig_range_ns=[],
        time_domain_ns=(-1_000, 5_000),
        min_sep_ns=6_000,
    )
    assert len(result_1w) == n_fp

    # two windows: 2 x n_fp RC events
    result_2w = spms_pars.get_rc_library(
        evt_file,
        lookup,
        ext_trig_range_ns=[(1_000, 7_000), (14_000, 20_000)],
        ge_trig_range_ns=[],
        time_domain_ns=(-1_000, 5_000),
        min_sep_ns=6_000,
    )
    assert len(result_2w) == 2 * n_fp

    # the rawid for the same source event is identical across windows
    # (window 0: indices 0..n_fp-1, window 1: indices n_fp..2*n_fp-1)
    if n_fp > 0:
        assert ak.to_list(result_2w.rawid[:n_fp]) == ak.to_list(result_2w.rawid[n_fp:])


def test_get_chunk_rc_data_returns_exact_size(legend_testdata):
    evt_file = legend_testdata[EVT_FILE]
    rc_evt_files = [str(evt_file)]
    lookup = spms_pars.build_rc_evt_index_lookup(rc_evt_files)

    chunk = spms_pars.get_chunk_rc_data(rc_evt_files, {}, 10, lookup)
    assert len(chunk) == 10


def test_get_chunk_rc_data_carryover(legend_testdata):
    """Leftover events from a file are reused in the next chunk."""
    evt_file = legend_testdata[EVT_FILE]
    rc_evt_files = [str(evt_file)]
    lookup = spms_pars.build_rc_evt_index_lookup(rc_evt_files)

    state: dict = {}
    chunk1 = spms_pars.get_chunk_rc_data(rc_evt_files, state, 10, lookup)
    assert len(chunk1) == 10
    # if the file had more events than 10, carryover should be populated
    lib_size = len(spms_pars.get_rc_library(evt_file, lookup))
    if lib_size > 10:
        assert state.get("carryover") is not None

    chunk2 = spms_pars.get_chunk_rc_data(rc_evt_files, state, 10, lookup)
    assert len(chunk2) == 10
