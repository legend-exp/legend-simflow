from __future__ import annotations

import awkward as ak
import lh5
import numpy as np
import pyg4ometry
import reboost.hpge

from legendsimflow import reboost as rutils
from legendsimflow import spms_pars


def test_psd_stuff(legend_testdata):
    dt_map = reboost.hpge.load_hpge_drift_time_maps(
        legend_testdata["lh5/V00048A-drift-time-maps-xtal-axes.lh5"],
        "V00048A",
        bounds_error=False,
    )

    xloc = [
        [0.162, 0.162, 0.162, 0.162, 0.162, 0.162, 0.162, 0.162, 0.162],
        [0.183, 0.183],
        [0.197, 0.197, 0.197, 0.197, 0.197, 0.197, 0.197, 0.197],
        [0.201, 0.201, 0.201, 0.201, 0.201, 0.201, 0.201, 0.201],
        [0.213, 0.213, 0.213, 0.213, 0.213],
        [0.176, 0.176, 0.176, 0.176, 0.176, 0.176, 0.177],
        [0.193, 0.193, 0.193, 0.193, 0.193, 0.193],
        [0.201, 0.201, 0.201],
        [0.157, 0.157, 0.157, 0.157, 0.157],
        [0.164, 0.164, 0.164],
    ]

    yloc = [
        [0.107, 0.107, 0.107, 0.107, 0.107, 0.107, 0.107, 0.107, 0.107],
        [0.0906, 0.0906],
        [0.122, 0.122, 0.122, 0.122, 0.122, 0.122, 0.122, 0.122],
        [0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11],
        [0.131, 0.131, 0.131, 0.131, 0.131],
        [0.0946, 0.0946, 0.0946, 0.0947, 0.0947, 0.0947, 0.0947],
        [0.157, 0.157, 0.157, 0.157, 0.157, 0.157],
        [0.145, 0.145, 0.145],
        [0.124, 0.124, 0.124, 0.124, 0.124],
        [0.144, 0.144, 0.144],
    ]

    zloc = [
        [0.535, 0.535, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536],
        [0.509, 0.509],
        [0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538],
        [0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533],
        [0.52, 0.52, 0.52, 0.52, 0.52],
        [0.506, 0.506, 0.506, 0.507, 0.507, 0.507, 0.507],
        [0.557, 0.557, 0.557, 0.557, 0.557, 0.557],
        [0.516, 0.516, 0.516],
        [0.519, 0.519, 0.519, 0.519, 0.519],
        [0.54, 0.54, 0.54],
    ]

    edep = ak.Array(
        [
            [130, 97.2, 233, 179, 133, 85.9, 173, 129, 331],
            [342, 125],
            [42.7, 75.9, 92.6, 261, 193, 49.6, 118, 249],
            [179, 35.6, 72.4, 174, 167, 241, 134, 126],
            [173, 37.3, 184, 156, 259],
            [13.9, 29.9, 91.1, 107, 106, 193, 70.8],
            [218, 13.2, 111, 111, 143, 118],
            [84.3, 56.9, 314],
            [63.7, 17.8, 92.4, 213, 186],
            [51, 95.3, 312],
        ]
    )

    chunk = ak.Array({"xloc": xloc, "yloc": yloc, "zloc": zloc})
    det_loc = pyg4ometry.gdml.Defines.Position(
        "Position",
        183.415,
        125.070,
        490.044,
        unit="mm",
    )

    dt = reboost.hpge.drift_time_crystal_axes(
        chunk.xloc,
        chunk.yloc,
        chunk.zloc,
        dt_map,
        coord_offset=det_loc,
    )

    assert ak.all((dt > 0) & (dt < 3000))

    pars = {
        "amax": 852.0,
        "mu": 0.0999,
        "sigma": 54.5,
        "tail_fraction": 0.223,
        "tau": 507.0,
        "high_tail_fraction": 0.00611,
        "high_tau": 593.0,
    }
    amax = rutils.hpge_max_current(edep, dt, pars)

    assert ak.all((amax > 0) & (amax < 3000))


def test_process_spms_windows_basic():
    """Test _process_spms_windows with simple data."""
    # Create mock SiPM data
    spms = ak.Array(
        {
            "t0": [
                [2000, 3000, 7000, 8000, 9000]
            ],  # Two hits in first window, two in second
            "energy": [[1.0, 2.0, 5.0, 3.0, 4.0]],
        }
    )

    win_ranges = [(1000, 14000)]  # Single range covering both windows
    time_domain_ns = (-1000, 5000)  # 6000 ns window
    min_sep_ns = 1000

    time = ak.flatten(spms.t0, axis=-1)
    energy = ak.flatten(spms.energy, axis=-1)
    npe, t0 = spms_pars._process_spms_windows(
        time, energy, win_ranges, time_domain_ns, min_sep_ns
    )

    # Should extract 2 windows: first (1000-7000), second (8000-14000)
    # Hits at 2000, 3000 should be included in the first window
    # Hit at 7000 is not included in the first window because of half-open interval convention: (spms.t0 >= wstart) & (spms.t0 < wend)
    # Hits at 8000, 9000 should be included in the second window
    # (lower edge is inclusive: t0 >= wstart)
    assert len(npe) == 4
    assert ak.sum(npe) == 10.0

    # Check that times are shifted to time_domain_ns
    flat_t0 = t0
    assert ak.all(flat_t0 >= time_domain_ns[0])
    assert ak.all(flat_t0 <= time_domain_ns[1])


def _get_rc_library(evt_file: str, **kwargs) -> ak.Array:
    lookup = spms_pars.build_rc_evt_index_lookup([evt_file])
    return spms_pars.get_rc_library(evt_file, lookup, **kwargs)


def test_forced_trigger_library_basic(legend_testdata):
    """Test get_random_coincidences_library with real event data."""
    f_evt = legend_testdata["lh5/l200-p16-r008-ssc-20251006T205904Z-tier_evt.lh5"]

    # Test with default parameters
    result = _get_rc_library(f_evt)

    # Check returned structure
    assert "rawid" in result.fields
    assert "npe" in result.fields
    assert "t0" in result.fields

    # Check that we got some data
    assert len(result) > 0

    # Check that times are in expected domain
    flat_t0 = ak.flatten(result.t0)
    if len(flat_t0) > 0:
        assert ak.all(flat_t0 >= -1000)  # Default time_domain_ns[0]
        assert ak.all(flat_t0 <= 5000)  # Default time_domain_ns[1]

    # Check that npe values are non-negative
    flat_npe = ak.flatten(result.npe)
    assert ak.all(flat_npe >= 0)


def test_forced_trigger_library_custom_time_domain(legend_testdata):
    """Test get_random_coincidences_library with custom time domain."""
    f_evt = legend_testdata["lh5/l200-p16-r008-ssc-20251006T205904Z-tier_evt.lh5"]

    # Use custom time domain
    result = _get_rc_library(f_evt, time_domain_ns=(-500, 3000))

    # Check that times are in custom domain
    flat_t0 = ak.flatten(result.t0)
    if len(flat_t0) > 0:
        assert ak.all(flat_t0 >= -500)
        assert ak.all(flat_t0 <= 3000)


def test_forced_trigger_library_custom_ranges(legend_testdata):
    """Test get_random_coincidences_library with custom window ranges."""
    f_evt = legend_testdata["lh5/l200-p16-r008-ssc-20251006T205904Z-tier_evt.lh5"]

    # Use custom ranges
    result = _get_rc_library(
        f_evt,
        ext_trig_range_ns=[(1000, 2000)],  # Single smaller range
        ge_trig_range_ns=[(1000, 2000)],  # Single smaller range
    )

    assert "npe" in result.fields
    assert "t0" in result.fields
    if len(result) > 0:
        flat_t0 = ak.flatten(result.t0)
        assert ak.all(flat_t0 >= -1000)
        assert ak.all(flat_t0 <= 5000)


def test_forced_trigger_library_rawid_consistency(legend_testdata):
    """Test reproducibility for repeated calls with identical inputs."""
    f_evt = legend_testdata["lh5/l200-p16-r008-ssc-20251006T205904Z-tier_evt.lh5"]

    result_1 = _get_rc_library(f_evt)
    result_2 = _get_rc_library(f_evt)

    assert len(result_1) == len(result_2)
    assert ak.to_list(result_1.rawid) == ak.to_list(result_2.rawid)
    assert ak.to_list(result_1.npe) == ak.to_list(result_2.npe)
    assert ak.to_list(result_1.t0) == ak.to_list(result_2.t0)


def test_forced_trigger_library_structure(legend_testdata):
    """Test the structure of returned data from get_random_coincidences_library."""
    f_evt = legend_testdata["lh5/l200-p16-r008-ssc-20251006T205904Z-tier_evt.lh5"]

    result = _get_rc_library(f_evt)

    # rawid, npe and t0 should have matching outer lengths
    assert len(result.rawid) == len(result.npe) == len(result.t0)


def test_forced_trigger_library_evt_number(legend_testdata):
    """Test the structure of returned data from get_random_coincidences_library."""
    f_evt = legend_testdata["lh5/l200-p16-r008-ssc-20251006T205904Z-tier_evt.lh5"]

    result = _get_rc_library(f_evt)

    evt = lh5.read_as("evt", f_evt, "ak")

    is_forced = evt.trigger.is_forced
    is_geds_trig = evt.coincident.geds
    is_muon = evt.coincident.muon_offline
    is_pulser = evt.coincident.puls  # codespell:ignore puls

    mask_forced_pulser = (is_forced | is_pulser) & ~is_muon
    num_forced_pulser = len(evt[mask_forced_pulser])

    mask_geds = is_geds_trig & ~is_muon
    num_geds = len(evt[mask_geds])

    # output contains only hits in requested windows from trigger-selected events
    total_triggers = num_forced_pulser + num_geds
    # If there are no eligible triggers, the RC library should be empty;
    # if there are eligible triggers, we expect at least one extracted window.
    if total_triggers == 0:
        assert len(result) == 0
    else:
        assert len(result) > 0
        # When entries are present, npe and t0 should have matching outer lengths
        assert len(result) == len(result.rawid) == len(result.npe) == len(result.t0)


def test_forced_trigger_library_num_processed_files(legend_testdata):
    f_evt = legend_testdata["lh5/l200-p16-r008-ssc-20251006T205904Z-tier_evt.lh5"]

    r1 = _get_rc_library(f_evt)
    r2 = _get_rc_library(f_evt)
    r3 = ak.concatenate(
        [
            _get_rc_library(f_evt),
            _get_rc_library(f_evt),
        ]
    )

    assert len(r1) == len(r2)
    assert len(r3) == 2 * len(r1)


def test_gauss_smear_output_type_and_shape():
    """Test gauss_smear returns an ak.Array with the same length as input."""
    arr_true = ak.Array([1.0, 2.0, 3.0, 4.0, 5.0])
    arr_reso = ak.Array([0.1, 0.1, 0.1, 0.1, 0.1])

    result = rutils.gauss_smear(arr_true, arr_reso)

    assert isinstance(result, ak.Array)
    assert len(result) == len(arr_true)


def test_gauss_smear_non_negative_for_positive_input():
    """Test gauss_smear does not return negative values for positive inputs."""
    arr_true = ak.Array(np.ones(1000))
    arr_reso = ak.Array(np.full(1000, 0.5))

    result = rutils.gauss_smear(arr_true, arr_reso)

    assert ak.all(result >= 0)


def test_gauss_smear_zero_replaced_by_tiny():
    """Test gauss_smear replaces near-zero smeared results with a tiny positive value."""
    # Input of zero with small resolution: smeared value may go negative and
    # should be replaced by np.finfo(float).tiny
    arr_true = ak.Array([0.0])
    arr_reso = ak.Array([0.0])

    result = rutils.gauss_smear(arr_true, arr_reso)

    assert isinstance(result, ak.Array)
    assert len(result) == 1
    assert result[0] >= 0
