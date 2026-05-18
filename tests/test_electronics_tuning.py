from __future__ import annotations

import numpy as np
import pytest

from legendsimflow.electronics_tuning import (
    build_cost_function,
    compute_rms_in_slice,
    select_ideal_wfs_in_slice,
)
from legendsimflow.superpulses import Slice, Superpulse


def _make_step_waveforms(dt, drift_times):
    n_samples = int(max(drift_times) / dt) + 200
    wfs = np.zeros((len(drift_times), n_samples))
    for i, t in enumerate(drift_times):
        idx = int(t / dt)
        wfs[i, idx:] = 1.0
    return wfs


def _make_superpulse(sl, current_wf, current_time_axis):
    """Minimal Superpulse with only the fields used by compute_rms_in_slice."""
    n = len(current_wf)
    return Superpulse(
        charge_wf=np.zeros(n),
        current_wf=current_wf,
        charge_time_axis=current_time_axis,
        current_time_axis=current_time_axis,
        slice=sl,
        detector="test",
        n_events_preliminary=100,
        n_events_final=100,
    )


def test_select_ideal_wfs_in_slice():
    dt = 1.0
    wfs = _make_step_waveforms(dt, [500, 1000, 1500, 2000])
    sl = Slice(energy_range=(0, 1e6), drift_time_range=(800, 1600))

    selected = select_ideal_wfs_in_slice(wfs, dt, sl)
    assert selected.shape[0] == 2


def test_select_ideal_wfs_in_slice_empty():
    dt = 1.0
    wfs = _make_step_waveforms(dt, [500])
    sl = Slice(energy_range=(0, 1e6), drift_time_range=(2000, 3000))

    selected = select_ideal_wfs_in_slice(wfs, dt, sl)
    assert len(selected) == 0


def test_select_ideal_wfs_in_slice_nan():
    dt = 1.0
    wfs = np.full((3, 200), np.nan)
    sl = Slice(energy_range=(0, 1e6), drift_time_range=(0, 1000))

    selected = select_ideal_wfs_in_slice(wfs, dt, sl)
    assert len(selected) == 0


def test_compute_rms_identical():
    sl = Slice(energy_range=(0, 1e6), drift_time_range=(900, 1100))
    wf = np.sin(np.linspace(0, 2 * np.pi, 200))
    time = np.arange(200, dtype=float)
    sp = _make_superpulse(sl, wf, time)

    rms = compute_rms_in_slice(wf, time, sp)
    assert rms == pytest.approx(0.0, abs=1e-12)


def test_compute_rms_with_offset():
    sl = Slice(energy_range=(0, 1e6), drift_time_range=(900, 1100))
    time = np.arange(100, dtype=float)
    data_wf = np.ones(100)
    sim_wf = data_wf + 0.5
    sp = _make_superpulse(sl, data_wf, time)

    rms = compute_rms_in_slice(sim_wf, time, sp)
    assert rms == pytest.approx(0.5, abs=1e-12)


def test_compute_rms_comparison_window():
    sl = Slice(energy_range=(0, 1e6), drift_time_range=(900, 1100))
    time = np.arange(100, dtype=float)
    data_wf = np.zeros(100)
    sim_wf = np.zeros(100)
    # Offset only outside the window
    sim_wf[:20] = 1.0
    sp = _make_superpulse(sl, data_wf, time)

    rms_full = compute_rms_in_slice(sim_wf, time, sp)
    rms_windowed = compute_rms_in_slice(sim_wf, time, sp, comparison_window=(30, 90))

    assert rms_full > 0
    assert rms_windowed == pytest.approx(0.0, abs=1e-12)


def test_compute_rms_no_overlap_raises():
    sl = Slice(energy_range=(0, 1e6), drift_time_range=(900, 1100))
    sp = _make_superpulse(sl, np.ones(10), np.arange(10, dtype=float))

    # Simulation on a completely disjoint time range
    with pytest.raises(ValueError, match="no valid samples"):
        compute_rms_in_slice(np.ones(10), np.arange(100, 110, dtype=float), sp)


def test_build_cost_function_penalty():
    """Non-positive parameters must return the penalty value."""
    sl = Slice(energy_range=(0, 1e6), drift_time_range=(900, 1100))
    sp = _make_superpulse(sl, np.zeros(10), np.arange(10, dtype=float))
    dummy_wfs = np.zeros((2, 100))
 
    cost = build_cost_function(
        {sl: dummy_wfs}, {sl: sp}, dt=1.0, alignment_idx=50, nsamples_output=10,
    )
    assert cost(-1, 30) == 1e6
    assert cost(5, -1) == 1e6
    assert cost(0, 30) == 1e6


@pytest.fixture
def cost_fixture():
    """Cost function and true parameters built from synthetic step waveforms."""
    from legendsimflow import psl
    sl = Slice(energy_range=(0, 1e6), drift_time_range=(900, 1100))
    dt = 1.0
    n_samples = 6000
    alignment_idx = n_samples // 2
    n_out = 3000

    ideal_wfs = np.zeros((5, n_samples))
    for i, t_step in enumerate([960, 980, 1000, 1020, 1040]):
        ideal_wfs[i, t_step:] = 1.0

    sigma_true, tau_true = 5.0, 40.0
    rf_true = psl.build_electronics_response_kernel(
        dt, mu_bandwidth=0.0, sigma_bandwidth=sigma_true, tau_rc=tau_true,
    )
    processed, _ = psl.process_ideal_waveforms(
        ideal_wfs, rf_true, dt, alignment_idx, n_out,
        mw_pars=psl.MW_PARS, dt_data=psl.DT_DATA,
    )
    data_wf = np.mean(processed, axis=0)
    data_time = (np.arange(n_out) - alignment_idx) * dt
    sp = _make_superpulse(sl, data_wf, data_time)

    cost = build_cost_function(
        {sl: ideal_wfs}, {sl: sp}, dt, alignment_idx, n_out,
    )
    return cost, sigma_true, tau_true


def test_build_cost_function_at_truth(cost_fixture):
    """Cost is near zero when evaluated at the true parameters."""
    cost, sigma_true, tau_true = cost_fixture
    assert cost(sigma_true, tau_true) < 1e-6


def test_build_cost_function_off_truth(cost_fixture):
    """Cost increases away from the true parameters."""
    cost, sigma_true, tau_true = cost_fixture
    assert cost(20.0, 100.0) > cost(sigma_true, tau_true)