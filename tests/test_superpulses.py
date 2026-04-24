from __future__ import annotations

import awkward as ak
import numpy as np
import pytest
from lgdo import Array, Scalar, Struct, lh5

from legendsimflow.superpulses import (
    Slice,
    Superpulse,
    apply_chi2_cut,
    compute_chi2_vs_superpulse,
    compute_superpulse,
    select_data_in_slice,
    select_detector_events,
    write_superpulses_to_lh5,
)

# ---------------------------------------------------------------------------
# Helpers / shared fixtures
# ---------------------------------------------------------------------------


def _make_superpulse(
    n_charge=100,
    n_current=80,
    energy_range=(1500.0, 2000.0),
    drift_time_range=(900.0, 1100.0),
    detector="V03422A",
    n_preliminary=50,
    n_final=40,
):
    """Return a minimal but valid Superpulse for reuse across tests."""
    rng = np.random.default_rng(0)
    return Superpulse(
        charge_wf=rng.random(n_charge),
        current_wf=rng.random(n_current),
        charge_time_axis=np.linspace(-500.0, 500.0, n_charge),
        current_time_axis=np.linspace(-500.0, 500.0, n_current),
        slice=Slice(energy_range, drift_time_range),
        detector=detector,
        n_events_preliminary=n_preliminary,
        n_events_final=n_final,
    )


def _make_evt_data(
    n_events=20,
    detector_names=None,
    hit_idxs=None,
    multiplicity=1,
    forced=False,
    puls=False,
    muon=False,
    muon_offline=False,
    is_good_channel=True,
    is_bb_like=True,
    spms_energy_sum=100.0,
    low_aoe=0.0,
    high_aoe=0.0,
    energy=1750.0,
    first_t0=200.0,
):
    """
    Build a minimal mock EVT awkward array.

    All per-detector arrays are wrapped in a length-1 inner dimension to
    reproduce the multiplicity==1 guarantee of the real EVT tier.
    """
    if detector_names is None:
        detector_names = [["V03422A"]] * n_events
    if hit_idxs is None:
        hit_idxs = [[i] for i in range(n_events)]

    return ak.Array(
        {
            "trigger": {"is_forced": [forced] * n_events},
            "coincident": {
                "puls": [puls] * n_events,
                "muon": [muon] * n_events,
                "muon_offline": [muon_offline] * n_events,
            },
            "geds": {
                "detector_name": detector_names,
                "hit_idx": hit_idxs,
                "multiplicity": [multiplicity] * n_events,
                "quality": {
                    "is_good_channel": [[is_good_channel]] * n_events,
                    "is_bb_like": [is_bb_like] * n_events,
                },
                "psd": {
                    "low_aoe": {"value": [[low_aoe]] * n_events},
                    "high_aoe": {"value": [[high_aoe]] * n_events},
                },
                "energy": [[energy]] * n_events,
            },
            "spms": {
                "energy_sum": [spms_energy_sum] * n_events,
                "first_t0": [first_t0] * n_events,
            },
        }
    )


# ===========================================================================
# Slice
# ===========================================================================


def test_slice_properties():
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    assert sl.energy_center == pytest.approx(1750.0)
    assert sl.drift_time_center == pytest.approx(1000.0)


def test_slice_str():
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    s = str(sl)
    assert "1500" in s
    assert "2000" in s
    assert "900" in s
    assert "1100" in s


def test_slice_is_hashable():
    sl1 = Slice((1500.0, 2000.0), (900.0, 1100.0))
    sl2 = Slice((1500.0, 2000.0), (900.0, 1100.0))
    d = {sl1: "value"}
    # same content -> same hash -> key lookup works
    assert d[sl2] == "value"


def test_slice_is_frozen():
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    with pytest.raises((AttributeError, TypeError)):
        sl.energy_range = (0.0, 1.0)


# ===========================================================================
# Superpulse construction and validation
# ===========================================================================


def test_superpulse_init_valid():
    sp = _make_superpulse()
    assert sp.detector == "V03422A"
    assert sp.n_events_preliminary == 50
    assert sp.n_events_final == 40


def test_superpulse_mismatched_charge_length():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="charge_wf"):
        Superpulse(
            charge_wf=rng.random(100),
            current_wf=rng.random(80),
            charge_time_axis=np.linspace(-500, 500, 99),  # wrong length
            current_time_axis=np.linspace(-500, 500, 80),
            slice=Slice((1500.0, 2000.0), (900.0, 1100.0)),
            detector="V03422A",
            n_events_preliminary=50,
            n_events_final=40,
        )


def test_superpulse_mismatched_current_length():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="current_wf"):
        Superpulse(
            charge_wf=rng.random(100),
            current_wf=rng.random(80),
            charge_time_axis=np.linspace(-500, 500, 100),
            current_time_axis=np.linspace(-500, 500, 50),  # wrong length
            slice=Slice((1500.0, 2000.0), (900.0, 1100.0)),
            detector="V03422A",
            n_events_preliminary=50,
            n_events_final=40,
        )


def test_superpulse_final_exceeds_preliminary():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="n_events_final"):
        Superpulse(
            charge_wf=rng.random(100),
            current_wf=rng.random(80),
            charge_time_axis=np.linspace(-500, 500, 100),
            current_time_axis=np.linspace(-500, 500, 80),
            slice=Slice((1500.0, 2000.0), (900.0, 1100.0)),
            detector="V03422A",
            n_events_preliminary=10,
            n_events_final=20,  # more than preliminary
        )


def test_superpulse_repr():
    sp = _make_superpulse()
    r = repr(sp)
    assert "V03422A" in r
    assert "40/50" in r


# ===========================================================================
# Superpulse.to_lgdo
# ===========================================================================


def test_to_lgdo_returns_struct():
    sp = _make_superpulse()
    result = sp.to_lgdo()
    assert isinstance(result, Struct)


def test_to_lgdo_fields_present():
    sp = _make_superpulse()
    result = sp.to_lgdo()
    expected_keys = {
        "charge_wf",
        "current_wf",
        "charge_time_axis",
        "current_time_axis",
        "dt_center",
        "dt_lo",
        "dt_hi",
        "e_lo",
        "e_hi",
        "detector",
        "n_events_preliminary",
        "n_events_final",
    }
    assert expected_keys == set(result.keys())


def test_to_lgdo_array_types():
    sp = _make_superpulse()
    result = sp.to_lgdo()
    assert isinstance(result["charge_wf"], Array)
    assert isinstance(result["current_wf"], Array)
    assert isinstance(result["charge_time_axis"], Array)
    assert isinstance(result["current_time_axis"], Array)


def test_to_lgdo_scalar_types():
    sp = _make_superpulse()
    result = sp.to_lgdo()
    for key in (
        "dt_center",
        "dt_lo",
        "dt_hi",
        "e_lo",
        "e_hi",
        "detector",
        "n_events_preliminary",
        "n_events_final",
    ):
        assert isinstance(result[key], Scalar), f"{key} should be Scalar"


def test_to_lgdo_scalar_values():
    sp = _make_superpulse(
        energy_range=(1500.0, 2000.0),
        drift_time_range=(900.0, 1100.0),
        n_preliminary=50,
        n_final=40,
    )
    result = sp.to_lgdo()
    assert result["dt_center"].value == pytest.approx(1000.0)
    assert result["dt_lo"].value == pytest.approx(900.0)
    assert result["dt_hi"].value == pytest.approx(1100.0)
    assert result["e_lo"].value == pytest.approx(1500.0)
    assert result["e_hi"].value == pytest.approx(2000.0)
    assert result["n_events_preliminary"].value == 50
    assert result["n_events_final"].value == 40
    assert result["detector"].value == "V03422A"


def test_to_lgdo_time_axis_units():
    sp = _make_superpulse()
    result = sp.to_lgdo()
    assert result["charge_time_axis"].attrs.get("units") == "ns"
    assert result["current_time_axis"].attrs.get("units") == "ns"


def test_to_lgdo_waveform_values():
    sp = _make_superpulse(n_charge=100, n_current=80)
    result = sp.to_lgdo()
    np.testing.assert_array_equal(result["charge_wf"].view_as("np"), sp.charge_wf)
    np.testing.assert_array_equal(result["current_wf"].view_as("np"), sp.current_wf)


# ===========================================================================
# select_detector_events
# ===========================================================================


def test_select_detector_events_basic():
    evt = _make_evt_data(n_events=10, detector_names=[["V03422A"]] * 10)
    result = select_detector_events(evt, "V03422A")
    assert len(result) == 10


def test_select_detector_events_filters_other_detector():
    # Mix 6 events for V03422A and 4 for B00035B
    detector_names = [["V03422A"]] * 6 + [["B00035B"]] * 4
    evt = _make_evt_data(n_events=10, detector_names=detector_names)
    result = select_detector_events(evt, "V03422A")
    assert len(result) == 6


def test_select_detector_events_none_match():
    evt = _make_evt_data(n_events=5, detector_names=[["B00035B"]] * 5)
    result = select_detector_events(evt, "V03422A")
    assert len(result) == 0


# ===========================================================================
# select_data_in_slice
# ===========================================================================


def test_select_data_in_slice_all_pass():
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    # energy=1750, drift_time must be added manually; we add it as a field
    evt = _make_evt_data(n_events=10, energy=1750.0)
    evt = ak.with_field(evt, [[1000.0]] * 10, "drift_time")
    result = select_data_in_slice(evt, sl)
    assert len(result) == 10


def test_select_data_in_slice_energy_out_of_range():
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    evt = _make_evt_data(n_events=5, energy=2500.0)  # above e_hi
    evt = ak.with_field(evt, [[1000.0]] * 5, "drift_time")
    result = select_data_in_slice(evt, sl)
    assert len(result) == 0


def test_select_data_in_slice_drift_time_out_of_range():
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    evt = _make_evt_data(n_events=5, energy=1750.0)
    evt = ak.with_field(evt, [[500.0]] * 5, "drift_time")  # below dt_lo
    result = select_data_in_slice(evt, sl)
    assert len(result) == 0


def test_select_data_in_slice_boundary_inclusive():
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    # events exactly on the boundary edges should be included
    evt = _make_evt_data(n_events=2, energy=1500.0)
    evt = ak.with_field(evt, [[900.0], [1100.0]], "drift_time")
    result = select_data_in_slice(evt, sl)
    assert len(result) == 2


def test_select_data_in_slice_mixed():
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    # Build four events with different energy/drift_time combinations:
    # event 0: E=1750 dt=1000 -> pass
    # event 1: E=1750 dt=500  -> fail (dt out of range)
    # event 2: E=2500 dt=1000 -> fail (energy out of range)
    # event 3: E=1750 dt=1000 -> pass
    evt0 = _make_evt_data(n_events=1, energy=1750.0)
    evt0 = ak.with_field(evt0, [[1000.0]], "drift_time")
    evt1 = _make_evt_data(n_events=1, energy=1750.0)
    evt1 = ak.with_field(evt1, [[500.0]], "drift_time")
    evt2 = _make_evt_data(n_events=1, energy=2500.0)
    evt2 = ak.with_field(evt2, [[1000.0]], "drift_time")
    evt3 = _make_evt_data(n_events=1, energy=1750.0)
    evt3 = ak.with_field(evt3, [[1000.0]], "drift_time")
    evt = ak.concatenate([evt0, evt1, evt2, evt3])
    result = select_data_in_slice(evt, sl)
    # only events 0 and 3 pass
    assert len(result) == 2


# ===========================================================================
# compute_superpulse
# ===========================================================================


def test_compute_superpulse_returns_superpulse():
    rng = np.random.default_rng(0)
    n_events, n_charge, n_current = 20, 100, 80
    charge_times = np.linspace(-500.0, 500.0, n_charge)
    current_times = np.linspace(-500.0, 500.0, n_current)
    charge_wfs = rng.random((n_events, n_charge))
    current_wfs = rng.random((n_events, n_current))
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))

    sp = compute_superpulse(
        charge_times,
        current_times,
        charge_wfs,
        current_wfs,
        sl,
        "V03422A",
        n_events,
    )
    assert isinstance(sp, Superpulse)


def test_compute_superpulse_shape():
    rng = np.random.default_rng(0)
    n_events, n_charge, n_current = 20, 100, 80
    charge_times = np.linspace(-500.0, 500.0, n_charge)
    current_times = np.linspace(-500.0, 500.0, n_current)
    charge_wfs = rng.random((n_events, n_charge))
    current_wfs = rng.random((n_events, n_current))
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))

    sp = compute_superpulse(
        charge_times,
        current_times,
        charge_wfs,
        current_wfs,
        sl,
        "V03422A",
        n_events,
    )
    assert sp.charge_wf.shape == (n_charge,)
    assert sp.current_wf.shape == (n_current,)


def test_compute_superpulse_is_mean():
    # With all waveforms identical, the superpulse must equal that waveform
    n_events, n_samples = 10, 50
    template = np.linspace(0.0, 1.0, n_samples)
    charge_wfs = np.tile(template, (n_events, 1))
    current_wfs = np.tile(template, (n_events, 1))
    times = np.linspace(-500.0, 500.0, n_samples)
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))

    sp = compute_superpulse(
        times, times, charge_wfs, current_wfs, sl, "V03422A", n_events
    )
    np.testing.assert_allclose(sp.charge_wf, template)
    np.testing.assert_allclose(sp.current_wf, template)


def test_compute_superpulse_nanmean():
    # NaN in one waveform should not corrupt the average
    n_events, n_samples = 4, 50
    times = np.linspace(-500.0, 500.0, n_samples)
    charge_wfs = np.ones((n_events, n_samples))
    charge_wfs[0, 10] = np.nan  # inject one NaN
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))

    sp = compute_superpulse(
        times, times, charge_wfs, charge_wfs.copy(), sl, "V03422A", n_events
    )
    assert np.isfinite(sp.charge_wf[10])
    assert sp.charge_wf[10] == pytest.approx(1.0)


def test_compute_superpulse_metadata():
    rng = np.random.default_rng(0)
    n_events, n_samples = 15, 50
    times = np.linspace(-500.0, 500.0, n_samples)
    charge_wfs = rng.random((n_events, n_samples))
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))

    sp = compute_superpulse(
        times,
        times,
        charge_wfs,
        charge_wfs,
        sl,
        "V03422A",
        n_events_preliminary=20,  # preliminary count passed separately
    )
    assert sp.n_events_final == n_events  # shape[0] of the input array
    assert sp.n_events_preliminary == 20
    assert sp.detector == "V03422A"
    assert sp.slice is sl


# ===========================================================================
# compute_chi2_vs_superpulse
# ===========================================================================


def _make_chi2_inputs(n_events=10, n_samples=100, seed=0):
    """Shared setup for chi2 tests.

    The time axis must extend well before -drift_time_range[1] so the
    auto-derived baseline mask (time < -1100 ns) selects some samples
    in fallback mode.
    """
    rng = np.random.default_rng(seed)
    times = np.linspace(-2000.0, 2000.0, n_samples)
    charge_wfs = rng.random((n_events, n_samples))
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    sp = compute_superpulse(
        times, times, charge_wfs, charge_wfs.copy(), sl, "V03422A", n_events
    )
    return charge_wfs, sp, n_events, n_samples


def test_compute_chi2_returns_array():
    charge_wfs, sp, n_events, _ = _make_chi2_inputs()
    bl_std = np.ones(n_events) * 2.0
    cuspEmax = np.ones(n_events) * 1000.0
    chi2 = compute_chi2_vs_superpulse(charge_wfs, sp, bl_std=bl_std, cuspEmax=cuspEmax)
    assert isinstance(chi2, np.ndarray)
    assert chi2.shape == (n_events,)


def test_compute_chi2_non_negative():
    charge_wfs, sp, n_events, _ = _make_chi2_inputs()
    bl_std = np.ones(n_events) * 2.0
    cuspEmax = np.ones(n_events) * 1000.0
    chi2 = compute_chi2_vs_superpulse(charge_wfs, sp, bl_std=bl_std, cuspEmax=cuspEmax)
    assert np.all(chi2 >= 0)


def test_compute_chi2_fallback_empty_baseline_raises():
    # Time axis that never satisfies time < -drift_time_range[1] = -1100 ns,
    # so the auto-derived baseline mask is all-False.
    # The module must raise ValueError rather than crash with a numpy warning.
    n_events, n_samples = 10, 100
    times = np.linspace(-500.0, 500.0, n_samples)
    charge_wfs = np.random.default_rng(0).random((n_events, n_samples))
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    sp = compute_superpulse(
        times, times, charge_wfs, charge_wfs.copy(), sl, "V03422A", n_events
    )
    with pytest.raises(ValueError, match="no baseline samples"):
        compute_chi2_vs_superpulse(charge_wfs, sp)


def test_compute_chi2_dsp_noise_mode():
    charge_wfs, sp, n_events, _ = _make_chi2_inputs()
    bl_std = np.ones(n_events) * 2.0
    cuspEmax = np.ones(n_events) * 1000.0
    chi2 = compute_chi2_vs_superpulse(charge_wfs, sp, bl_std=bl_std, cuspEmax=cuspEmax)
    assert chi2.shape == (n_events,)
    assert np.all(chi2 >= 0)


def test_compute_chi2_dsp_vs_fallback_shape():
    # Both modes must return arrays of the same shape.
    # Pass an explicit baseline_region_mask so the fallback does not depend
    # on the auto-derived mask (which requires a specific time axis range).
    charge_wfs, sp, n_events, n_samples = _make_chi2_inputs()
    bl_std = np.ones(n_events) * 2.0
    cuspEmax = np.ones(n_events) * 1000.0
    baseline_mask = np.zeros(n_samples, dtype=bool)
    baseline_mask[:20] = True
    chi2_dsp = compute_chi2_vs_superpulse(
        charge_wfs, sp, bl_std=bl_std, cuspEmax=cuspEmax
    )
    chi2_fallback = compute_chi2_vs_superpulse(
        charge_wfs, sp, baseline_region_mask=baseline_mask
    )
    assert chi2_dsp.shape == chi2_fallback.shape


def test_compute_chi2_identical_waveforms_is_small():
    # Waveforms identical to the superpulse should have chi2 near zero
    n_events, n_samples = 10, 100
    times = np.linspace(-500.0, 500.0, n_samples)
    template = np.linspace(0.0, 1.0, n_samples)
    charge_wfs = np.tile(template, (n_events, 1))
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    sp = compute_superpulse(
        times, times, charge_wfs, charge_wfs.copy(), sl, "V03422A", n_events
    )
    # Use a realistic noise level so the chi2 is well-defined
    bl_std = np.ones(n_events) * 2.0
    cuspEmax = np.ones(n_events) * 1000.0
    chi2 = compute_chi2_vs_superpulse(charge_wfs, sp, bl_std=bl_std, cuspEmax=cuspEmax)
    assert np.all(chi2 < 1e-6)


def test_compute_chi2_outlier_is_larger():
    # An outlier waveform should have higher chi2 than the rest
    rng = np.random.default_rng(42)
    n_events, n_samples = 10, 100
    times = np.linspace(-500.0, 500.0, n_samples)
    charge_wfs = np.tile(np.linspace(0.0, 1.0, n_samples), (n_events, 1))
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    sp = compute_superpulse(
        times, times, charge_wfs, charge_wfs.copy(), sl, "V03422A", n_events
    )
    # Replace last waveform with pure noise (guaranteed outlier)
    outlier_idx = -1
    charge_wfs_with_outlier = charge_wfs.copy()
    charge_wfs_with_outlier[outlier_idx] = rng.random(n_samples) * 10
    bl_std = np.ones(n_events) * 2.0
    cuspEmax = np.ones(n_events) * 1000.0
    chi2 = compute_chi2_vs_superpulse(
        charge_wfs_with_outlier, sp, bl_std=bl_std, cuspEmax=cuspEmax
    )
    assert chi2[outlier_idx] > chi2[:-1].max()


def test_compute_chi2_custom_baseline_mask():
    charge_wfs, sp, n_events, n_samples = _make_chi2_inputs()
    mask = np.zeros(n_samples, dtype=bool)
    mask[:20] = True  # first 20 samples as baseline
    chi2 = compute_chi2_vs_superpulse(charge_wfs, sp, baseline_region_mask=mask)
    assert chi2.shape == (n_events,)
    assert np.all(chi2 >= 0)


# ===========================================================================
# apply_chi2_cut
# ===========================================================================


def test_apply_chi2_cut_returns_correct_shapes():
    rng = np.random.default_rng(0)
    n_events, n_charge, n_current = 10, 100, 80
    charge_wfs = rng.random((n_events, n_charge))
    current_wfs = rng.random((n_events, n_current))
    chi2 = np.array([0.5, 4.0, 1.0, 2.5, 10.0, 0.8, 3.1, 0.9, 1.5, 6.0])
    threshold = 3.0

    golden_charge, golden_current, golden_idx = apply_chi2_cut(
        charge_wfs, current_wfs, chi2, threshold=threshold
    )
    n_expected = int((chi2 < threshold).sum())
    assert golden_charge.shape == (n_expected, n_charge)
    assert golden_current.shape == (n_expected, n_current)
    assert golden_idx.shape == (n_expected,)


def test_apply_chi2_cut_waveforms_match_indices():
    rng = np.random.default_rng(0)
    n_events = 5
    chi2 = np.array([0.5, 4.0, 1.0, 2.5, 10.0])
    charge_wfs = rng.random((n_events, 50))
    current_wfs = rng.random((n_events, 50))

    golden_charge, golden_current, golden_idx = apply_chi2_cut(
        charge_wfs, current_wfs, chi2, threshold=3.0
    )
    # The returned arrays must be the rows selected by golden_idx
    np.testing.assert_array_equal(golden_charge, charge_wfs[golden_idx])
    np.testing.assert_array_equal(golden_current, current_wfs[golden_idx])


def test_apply_chi2_cut_all_pass():
    rng = np.random.default_rng(0)
    n_events = 5
    chi2 = np.zeros(n_events)  # all below any positive threshold
    charge_wfs = rng.random((n_events, 50))
    current_wfs = rng.random((n_events, 50))

    _golden_charge, _golden_current, golden_idx = apply_chi2_cut(
        charge_wfs, current_wfs, chi2, threshold=3.0
    )
    assert len(golden_idx) == n_events


def test_apply_chi2_cut_none_pass():
    rng = np.random.default_rng(0)
    n_events = 5
    chi2 = np.full(n_events, 100.0)  # all above threshold
    charge_wfs = rng.random((n_events, 50))
    current_wfs = rng.random((n_events, 50))

    golden_charge, golden_current, golden_idx = apply_chi2_cut(
        charge_wfs, current_wfs, chi2, threshold=3.0
    )
    assert len(golden_idx) == 0
    assert golden_charge.shape == (0, 50)
    assert golden_current.shape == (0, 50)


def test_apply_chi2_cut_custom_threshold():
    rng = np.random.default_rng(0)
    n_events = 5
    chi2 = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    charge_wfs = rng.random((n_events, 50))
    current_wfs = rng.random((n_events, 50))

    _, _, idx_tight = apply_chi2_cut(charge_wfs, current_wfs, chi2, threshold=2.5)
    _, _, idx_loose = apply_chi2_cut(charge_wfs, current_wfs, chi2, threshold=4.5)

    assert len(idx_tight) == 2  # chi2 < 2.5: events 0 and 1
    assert len(idx_loose) == 4  # chi2 < 4.5: events 0, 1, 2, 3


# ===========================================================================
# write_superpulses_to_lh5
# ===========================================================================


def test_write_superpulses_creates_file(tmp_path):
    sp = _make_superpulse()
    output_path = str(tmp_path / "test_superpulses.lh5")
    superpulses = {sp.slice: sp}
    write_superpulses_to_lh5(superpulses, output_path, detector="V03422A")
    assert (tmp_path / "test_superpulses.lh5").exists()


def test_write_superpulses_lh5_structure(tmp_path):
    sp = _make_superpulse(
        drift_time_range=(900.0, 1100.0),
    )
    output_path = str(tmp_path / "test_superpulses.lh5")
    write_superpulses_to_lh5({sp.slice: sp}, output_path, detector="V03422A")

    # The expected group path is V03422A/dt_900_1100_ns
    result = lh5.read("V03422A/dt_900_1100_ns", output_path)
    assert isinstance(result, Struct)


def test_write_superpulses_lh5_field_values(tmp_path):
    sp = _make_superpulse(
        n_charge=100,
        n_current=80,
        drift_time_range=(900.0, 1100.0),
        energy_range=(1500.0, 2000.0),
        n_preliminary=50,
        n_final=40,
    )
    output_path = str(tmp_path / "test_superpulses.lh5")
    write_superpulses_to_lh5({sp.slice: sp}, output_path, detector="V03422A")

    result = lh5.read("V03422A/dt_900_1100_ns", output_path)

    np.testing.assert_array_equal(result["charge_wf"].view_as("np"), sp.charge_wf)
    np.testing.assert_array_equal(result["current_wf"].view_as("np"), sp.current_wf)
    assert result["dt_center"].value == pytest.approx(1000.0)
    assert result["e_lo"].value == pytest.approx(1500.0)
    assert result["e_hi"].value == pytest.approx(2000.0)
    assert result["n_events_preliminary"].value == 50
    assert result["n_events_final"].value == 40


def test_write_superpulses_multiple_slices(tmp_path):
    sl1 = Slice((1500.0, 2000.0), (900.0, 1100.0))
    sl2 = Slice((1500.0, 2000.0), (1100.0, 1300.0))
    sp1 = _make_superpulse(drift_time_range=(900.0, 1100.0))
    sp2 = _make_superpulse(drift_time_range=(1100.0, 1300.0))

    output_path = str(tmp_path / "multi.lh5")
    write_superpulses_to_lh5({sl1: sp1, sl2: sp2}, output_path, detector="V03422A")

    # Both groups must be present
    result1 = lh5.read("V03422A/dt_900_1100_ns", output_path)
    result2 = lh5.read("V03422A/dt_1100_1300_ns", output_path)
    assert isinstance(result1, Struct)
    assert isinstance(result2, Struct)
