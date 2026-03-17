"""Tests for spms_pars module."""

from __future__ import annotations

import awkward as ak
import pytest

from legendsimflow import spms_pars


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
