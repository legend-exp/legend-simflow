"""Tests for spms_pars module."""

from __future__ import annotations

import awkward as ak
import numpy as np
import pytest

from legendsimflow import spms_pars


@pytest.fixture
def mock_forced_trig_library():
    """Create a mock forced trigger library for testing."""
    # Create structured array with npe, t0, and rawid fields
    # 10 events, 3 SiPM channels each
    npe_data = ak.Array(
        [
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
            [[2, 3, 4], [5, 6, 7], [8, 9, 10]],
            [[3, 4, 5], [6, 7, 8], [9, 10, 11]],
            [[4, 5, 6], [7, 8, 9], [10, 11, 12]],
            [[5, 6, 7], [8, 9, 10], [11, 12, 13]],
            [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
            [[2, 2, 2], [3, 3, 3], [4, 4, 4]],
            [[3, 3, 3], [4, 4, 4], [5, 5, 5]],
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
            [[10, 20, 30], [40, 50, 60], [70, 80, 90]],
        ]
    )

    t0_data = ak.Array(
        [
            [[100.0, 200.0, 300.0], [150.0, 250.0, 350.0], [120.0, 220.0, 320.0]],
            [[105.0, 205.0, 305.0], [155.0, 255.0, 355.0], [125.0, 225.0, 325.0]],
            [[110.0, 210.0, 310.0], [160.0, 260.0, 360.0], [130.0, 230.0, 330.0]],
            [[115.0, 215.0, 315.0], [165.0, 265.0, 365.0], [135.0, 235.0, 335.0]],
            [[120.0, 220.0, 320.0], [170.0, 270.0, 370.0], [140.0, 240.0, 340.0]],
            [[50.0, 150.0, 250.0], [100.0, 200.0, 300.0], [75.0, 175.0, 275.0]],
            [[55.0, 155.0, 255.0], [105.0, 205.0, 305.0], [80.0, 180.0, 280.0]],
            [[60.0, 160.0, 260.0], [110.0, 210.0, 310.0], [85.0, 185.0, 285.0]],
            [[65.0, 165.0, 265.0], [115.0, 215.0, 315.0], [90.0, 190.0, 290.0]],
            [[70.0, 170.0, 270.0], [120.0, 220.0, 320.0], [95.0, 195.0, 295.0]],
        ]
    )

    # SiPM channel UIDs: 101, 102, 103
    rawid_array = np.array([101, 102, 103], dtype=np.int32)

    # Create the library structure like get_forced_trigger_library returns
    return ak.Array(
        {
            "npe": npe_data,
            "t0": t0_data,
            "rawid": ak.Array([rawid_array for _ in range(len(npe_data))]),
        }
    )


def test_extract_single_channel(mock_forced_trig_library):
    """Test extracting data for a single SiPM channel."""
    npe, t0 = spms_pars.get_sipm_rc_data(mock_forced_trig_library, "sipm", 101)

    # Should have data for all 10 events
    assert len(npe) == 10
    assert len(t0) == 10


def test_extract_all_channels(mock_forced_trig_library):
    """Test extracting data for all channels combined."""
    npe, t0 = spms_pars.get_sipm_rc_data(mock_forced_trig_library, "all", 0)

    # Should have flattened data but 10 events
    assert len(npe) == 10
    assert len(t0) == 10


def test_extract_different_channels(mock_forced_trig_library):
    """Test extracting data from different channels."""
    npe1, _ = spms_pars.get_sipm_rc_data(mock_forced_trig_library, "sipm", 101)
    npe2, _ = spms_pars.get_sipm_rc_data(mock_forced_trig_library, "sipm", 102)
    npe3, _ = spms_pars.get_sipm_rc_data(mock_forced_trig_library, "sipm", 103)

    # All should have same length but different values
    assert len(npe1) == len(npe2) == len(npe3) == 10
    # Values should be different
    assert ak.to_list(npe1) != ak.to_list(npe2)
    assert ak.to_list(npe2) != ak.to_list(npe3)


def test_extract_missing_sipm_raises_error(mock_forced_trig_library):
    """Test that missing SiPM UID raises ValueError."""
    with pytest.raises(ValueError, match="SiPM UID 999 not found"):
        spms_pars.get_sipm_rc_data(mock_forced_trig_library, "sipm", 999)


def test_extract_preserves_array_structure(mock_forced_trig_library):
    """Test that extracted data preserves awkward array structure."""
    npe, t0 = spms_pars.get_sipm_rc_data(mock_forced_trig_library, "sipm", 101)

    # Should be awkward arrays
    assert isinstance(npe, ak.Array)
    assert isinstance(t0, ak.Array)

    # Should contain valid data
    npe_list = ak.to_list(npe)
    t0_list = ak.to_list(t0)

    for npe_val, t0_val in zip(npe_list, t0_list, strict=True):
        assert isinstance(npe_val, list)
        assert isinstance(t0_val, list)
        assert len(npe_val) > 0
        assert len(t0_val) > 0


def test_extract_all_combines_channels(mock_forced_trig_library):
    """Test that 'all' mode properly flattens all channels."""
    npe_all, _ = spms_pars.get_sipm_rc_data(mock_forced_trig_library, "all", 0)

    npe_101, _ = spms_pars.get_sipm_rc_data(mock_forced_trig_library, "sipm", 101)
    npe_102, _ = spms_pars.get_sipm_rc_data(mock_forced_trig_library, "sipm", 102)
    npe_103, _ = spms_pars.get_sipm_rc_data(mock_forced_trig_library, "sipm", 103)

    # All data should be combination of three channels
    assert len(npe_all) == len(npe_101) == len(npe_102) == len(npe_103)
