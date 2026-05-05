from __future__ import annotations

import awkward as ak
import lgdo
import numpy as np
import pyg4ometry
import pytest
import reboost
import reboost.math.stats
from lgdo import lh5

from legendsimflow import reboost as rutils
from legendsimflow import spms_pars


def test_remage_hit_range(legend_testdata):
    f_stp = legend_testdata["remage/th228-full-optional-v0_13.lh5"]

    tcm = lh5.read_as("tcm", f_stp, library="ak")

    for det, uid in zip(
        ["det1", "det2", "scint1", "scint2", "optdet1", "optdet2"],
        [11, 12, 1, 2, 101, 102],
        strict=True,
    ):
        n_rows = lh5.read_n_rows(f"stp/{det}", f_stp)

        tcm = lh5.read_as("tcm", f_stp, "ak")

        assert rutils.get_remage_hit_range(tcm, det, uid, [0, len(tcm) - 1]) == (
            0,
            n_rows,
        )

        # divide into groups
        groups = [[0, 10], [11, 40], [41, 101], [102, len(tcm) - 1]]
        n = 0

        for group in groups:
            nen = rutils.get_remage_hit_range(tcm, det, uid, group)[1]

            if nen is not None:
                n += nen

        assert n == n_rows


def test_psd_stuff(legend_testdata):
    dt_map = {}
    for angle in ("000", "045"):
        dt_map[angle] = reboost.hpge.utils.get_hpge_rz_field(
            legend_testdata["lh5/V00048A-drift-time-maps-xtal-axes.lh5"],
            "V00048A",
            f"drift_time_{angle}_deg",
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

    dt = rutils.hpge_corrected_drift_time(
        chunk,
        dt_map,
        det_loc,
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


def test_cluster_photoelectrons_does_not_cross_subarrays():
    """Test that clustering does not merge elements across subarray boundaries."""
    times = ak.Array([[[0.0, 0.6], [0.7, 0.9]]])
    amps = ak.Array([[[1.0, 2.0], [3.0, 4.0]]])

    t_out, a_out = rutils.cluster_photoelectrons(times, amps, thr=1.0)

    assert ak.to_list(t_out) == [[[0.0], [0.7]]]
    assert ak.to_list(a_out) == [[[3.0], [7.0]]]


def test_cluster_photoelectrons_enforces_max_span():
    """Test that clusters respect the maximum time span threshold."""
    times = ak.Array([[0.0, 0.6, 1.1, 1.4, 2.3]])
    amps = ak.Array([[1.0, 2.0, 3.0, 4.0, 5.0]])

    t_out, a_out = rutils.cluster_photoelectrons(times, amps, thr=1.0)

    assert ak.to_list(t_out) == [[0.0, 1.1, 2.3]]
    assert ak.to_list(a_out) == [[3.0, 7.0, 5.0]]


def test_cluster_photoelectrons_empty_and_boundary():
    """Test clustering with empty arrays and exact boundary conditions."""
    times = ak.Array([[], [0.0, 1.0, 1.0001]])
    amps = ak.Array([[], [1.0, 2.0, 3.0]])

    t_out, a_out = rutils.cluster_photoelectrons(times, amps, thr=1.0)

    # [0.0, 1.0] spans exactly 1.0 -> same cluster; 1.0001 starts new
    assert ak.to_list(t_out) == [[], [0.0, 1.0001]]
    assert ak.to_list(a_out) == [[], [3.0, 3.0]]


def test_cluster_photoelectrons_mismatched_shapes():
    """Test that mismatched array shapes raise ValueError."""
    # Different nesting depths
    times_1d = ak.Array([0.0, 1.0, 2.0])
    amps_2d = ak.Array([[1.0, 2.0, 3.0]])

    with pytest.raises(ValueError, match="nesting depth"):
        rutils.cluster_photoelectrons(times_1d, amps_2d, thr=1.0)

    # Same nesting but different list lengths
    times = ak.Array([[0.0, 1.0], [2.0]])
    amps = ak.Array([[1.0], [2.0, 3.0]])

    with pytest.raises(ValueError, match="mismatched list lengths"):
        rutils.cluster_photoelectrons(times, amps, thr=1.0)


def test_smear_photoelectrons_shape_preservation():
    """Test that smear_photoelectrons preserves input array shape for 1D ragged arrays."""
    # Test with 1D ragged arrays (the intended use case)
    array_1d = ak.Array([[1.0, 2.0, 3.0], [4.0], [5.0, 6.0]])
    array_empty = ak.Array([[], [1.0, 2.0], []])

    rng = np.random.default_rng(42)
    result_1d = rutils.smear_photoelectrons(array_1d, fwhm_in_pe=0.5, rng=rng)
    # Check that the structure is preserved (number of elements per sublist)
    assert ak.num(result_1d, axis=1).to_list() == ak.num(array_1d, axis=1).to_list()
    assert len(result_1d) == len(array_1d)

    rng = np.random.default_rng(42)
    result_empty = rutils.smear_photoelectrons(array_empty, fwhm_in_pe=0.5, rng=rng)
    assert (
        ak.num(result_empty, axis=1).to_list() == ak.num(array_empty, axis=1).to_list()
    )
    assert len(result_empty) == len(array_empty)


def test_smear_photoelectrons_non_negativity():
    """Test that smear_photoelectrons clamps negative values to zero."""
    # Use large FWHM to increase probability of negative samples
    # With loc=1 and sigma=2/2.35482≈0.85, ~12% of samples would be negative
    array = ak.Array([np.ones(10000)])
    rng = np.random.default_rng(42)
    result = rutils.smear_photoelectrons(array, fwhm_in_pe=2.0, rng=rng)

    # Verify no negative values
    flat_result = ak.flatten(result, axis=None)
    assert ak.all(flat_result >= 0)

    # Verify that clamping actually occurred (some zeros should exist)
    assert ak.sum(flat_result == 0) > 0


def test_smear_photoelectrons_statistical_properties():
    """Test that smear_photoelectrons produces correct statistical distribution."""
    # Generate large sample to test statistical properties
    n_samples = 100000
    array = ak.Array([np.ones(n_samples)])
    fwhm = 0.4
    expected_sigma = fwhm / 2.35482

    rng = np.random.default_rng(42)
    result = rutils.smear_photoelectrons(array, fwhm_in_pe=fwhm, rng=rng)

    flat_result = ak.flatten(result, axis=None)
    # Filter out clamped zeros for statistical analysis of untruncated distribution
    non_zero = flat_result[flat_result > 0]

    # Mean should be close to 1 (allowing for small statistical fluctuation)
    # With small FWHM, very few values get clamped, so mean ≈ 1
    mean = ak.mean(non_zero)
    assert 0.99 < mean < 1.01

    # Standard deviation should match expected value (with tolerance)
    # Excluding clamped values gives us the true Gaussian sigma
    std = ak.std(non_zero)
    assert 0.95 * expected_sigma < std < 1.05 * expected_sigma


def test_smear_photoelectrons_reproducibility():
    """Test that smear_photoelectrons is reproducible with same seed."""
    array = ak.Array([[1.0, 2.0, 3.0, 4.0]])
    fwhm = 0.5

    rng1 = np.random.default_rng(123)
    result1 = rutils.smear_photoelectrons(array, fwhm_in_pe=fwhm, rng=rng1)

    rng2 = np.random.default_rng(123)
    result2 = rutils.smear_photoelectrons(array, fwhm_in_pe=fwhm, rng=rng2)

    assert ak.all(result1 == result2)


def test_smear_photoelectrons_default_rng():
    """Test that smear_photoelectrons works with default RNG (no rng parameter)."""
    array = ak.Array([[1.0, 2.0, 3.0]])

    # Should not raise an error
    result = rutils.smear_photoelectrons(array, fwhm_in_pe=0.5)

    # Should still preserve shape and non-negativity
    assert ak.num(result, axis=1).to_list() == ak.num(array, axis=1).to_list()
    assert ak.all(ak.flatten(result, axis=None) >= 0)


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


def test_get_remage_detector_uids(legend_testdata):
    """Test that get_remage_detector_uids returns a dict mapping UIDs to detector names."""
    f_stp = legend_testdata["remage/th228-full-optional-v0_13.lh5"]

    result = rutils.get_remage_detector_uids(f_stp)

    assert isinstance(result, dict)
    # all keys should be ints, all values should be strings
    assert all(isinstance(k, int) for k in result)
    assert all(isinstance(v, str) for v in result.values())
    # verify the known mapping from the test file
    assert result == {
        1: "scint1",
        2: "scint2",
        11: "det1",
        12: "det2",
        101: "optdet1",
        102: "optdet2",
    }


def test_make_output_chunk_scalar_t0():
    """Test make_output_chunk when chunk already has scalar t0 and evtid arrays."""
    chunk = lgdo.Table(size=3)
    chunk.add_field(
        "t0", lgdo.Array(np.array([100.0, 200.0, 300.0]), attrs={"units": "ns"})
    )
    chunk.add_field("evtid", lgdo.Array(np.array([1, 2, 3])))

    result = rutils.make_output_chunk(chunk)

    assert isinstance(result, lgdo.Table)
    assert len(result) == 3
    assert "t0" in result
    assert "evtid" in result
    np.testing.assert_array_equal(result.t0.nda, [100.0, 200.0, 300.0])
    np.testing.assert_array_equal(result.evtid.nda, [1, 2, 3])


def test_make_output_chunk_vector_time():
    """Test make_output_chunk when chunk has vector time and evtid (takes firsts)."""
    chunk = lgdo.Table(size=3)
    chunk.add_field(
        "time",
        lgdo.VectorOfVectors(
            ak.Array([[100.0, 200.0], [300.0], [400.0, 500.0, 600.0]]),
            attrs={"units": "ns"},
        ),
    )
    chunk.add_field("evtid", lgdo.VectorOfVectors(ak.Array([[1, 2], [3], [4, 5, 6]])))

    result = rutils.make_output_chunk(chunk)

    assert isinstance(result, lgdo.Table)
    assert len(result) == 3
    assert "t0" in result
    assert "evtid" in result
    # t0 should be the first time value of each event
    np.testing.assert_array_equal(result.t0.nda, [100.0, 300.0, 400.0])
    # evtid should be the first evtid of each event
    np.testing.assert_array_equal(result.evtid.nda, [1, 3, 4])


def test_write_chunk(tmp_path):
    """Test write_chunk writes data and creates soft links for new detectors."""
    outfile = tmp_path / "output.lh5"
    chunk = lgdo.Table(size=3)
    chunk.add_field(
        "t0", lgdo.Array(np.array([100.0, 200.0, 300.0]), attrs={"units": "ns"})
    )
    chunk.add_field("evtid", lgdo.Array(np.array([1, 2, 3])))

    # First write: creates file and soft link
    rutils.write_chunk(chunk, "hit/det001", outfile, 1)
    assert "hit/det001" in lh5.ls(outfile, "hit/")
    assert "hit/__by_uid__" in lh5.ls(outfile, "hit/")
    assert "hit/__by_uid__/det001" in lh5.ls(outfile, "hit/__by_uid__/")

    # Second write: appends to existing detector table (no new link)
    rutils.write_chunk(chunk, "hit/det001", outfile, 1)
    result = lh5.read("hit/det001", outfile)
    assert len(result) == 6  # appended 3 more rows

    # Third write: adds a new detector with its own link
    rutils.write_chunk(chunk, "hit/det002", outfile, 2)
    assert "hit/__by_uid__/det002" in lh5.ls(outfile, "hit/__by_uid__/")


def test_listoffset_chain_1d():
    """Test _listoffset_chain with a 1D list-of-values array."""
    arr = ak.Array([[1.0, 2.0], [3.0]])
    layout = ak.to_layout(arr)

    offsets_chain, content = rutils._listoffset_chain(layout)

    assert len(offsets_chain) == 1
    np.testing.assert_array_equal(offsets_chain[0], [0, 2, 3])
    assert isinstance(content, ak.contents.NumpyArray)


def test_listoffset_chain_nested():
    """Test _listoffset_chain with a doubly-nested array (3D)."""
    arr = ak.Array([[[1.0, 2.0], [3.0]], [[4.0, 5.0, 6.0]]])
    layout = ak.to_layout(arr)

    offsets_chain, content = rutils._listoffset_chain(layout)

    assert len(offsets_chain) == 2
    assert isinstance(content, ak.contents.NumpyArray)


def test_listoffset_chain_non_numpy_content():
    """Test _listoffset_chain raises TypeError when content is not NumpyArray."""
    # A list of records does not end in a NumpyArray
    arr = ak.Array([{"x": 1}, {"x": 2}])
    layout = ak.to_layout(arr)

    with pytest.raises(TypeError, match="NumpyArray"):
        rutils._listoffset_chain(layout)


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
