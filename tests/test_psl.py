from __future__ import annotations

import awkward as ak
import numpy as np
import pytest
from lgdo import Array, Scalar

from legendsimflow import psl


def test_build_electronics_response_kernel():
    # test first with all options
    dt = 10.0
    mu_bandwidth = 0.0
    sigma_bandwidth = 50.0
    tau_rc = 100.0

    kernel = psl.build_electronics_response_kernel(
        dt, mu_bandwidth, sigma_bandwidth, tau_rc
    )
    assert isinstance(kernel, np.ndarray)
    assert kernel.size > 0
    assert np.all(kernel >= 0)
    assert np.isclose(np.sum(kernel), 1.0, atol=1e-2)

    # 'full' mode results in a length of (len(rf_preamp) + len(rf_digi) - 1)
    assert len(kernel) == 600 * 2 - 1  # default kernel length

    # check the errors raised
    with pytest.raises(ValueError):
        psl.build_electronics_response_kernel(-1, mu_bandwidth, sigma_bandwidth, tau_rc)

    with pytest.raises(ValueError):
        psl.build_electronics_response_kernel(dt, mu_bandwidth, -10, tau_rc)

    with pytest.raises(ValueError):
        psl.build_electronics_response_kernel(dt, mu_bandwidth, sigma_bandwidth, -10)

    # now check with gaussian only

    dt = 10.0
    mu_bandwidth = 0.0
    sigma_bandwidth = 50.0
    tau_rc = 100.0

    kernel = psl.build_electronics_response_kernel(
        dt, mu_bandwidth, sigma_bandwidth, tau_rc, gaussian_only=True
    )
    assert isinstance(kernel, np.ndarray)
    assert kernel.size > 0
    assert np.all(kernel >= 0)
    assert np.isclose(np.sum(kernel), 1.0, atol=1e-2)

    assert len(kernel) == 600  # default kernel length for gaussian only


def test_convolve_waveform_with_kernel():
    # create a simple test waveform (delta function)
    wf = np.zeros((10, 1000))
    wf[:, 500] = 1.0

    dt = 10.0
    mu_bandwidth = 0.0
    sigma_bandwidth = 50.0
    tau_rc = 100.0

    kernel = psl.build_electronics_response_kernel(
        dt, mu_bandwidth, sigma_bandwidth, tau_rc
    )

    # convolve
    convolved = psl.apply_electronics_response(wf, kernel)

    assert isinstance(convolved, ak.Array)

    # shape should not be altered due to the slicing
    assert convolved.to_numpy().shape == wf.shape

    # test with a large batch size
    convolved_1000 = psl.apply_electronics_response(wf, kernel, batch_size=1000)

    # test with a small
    convolved_2 = psl.apply_electronics_response(wf, kernel, batch_size=2)

    assert convolved_1000.to_numpy().shape == wf.shape
    assert convolved_2.to_numpy().shape == wf.shape


def test_shift_array():
    # create a simple test waveform (delta function)
    wf = np.zeros((10, 1000))
    wf[:, 500] = 1.0

    shifted, peak_indices = psl.align_waveforms_to_peak(
        wf, alignment_idx=510, nsamples_output_current_wfs=1000
    )

    assert isinstance(shifted, np.ndarray)
    assert isinstance(peak_indices, np.ndarray)

    assert shifted.shape == wf.shape

    assert len(peak_indices) == wf.shape[0]


def test_align_waveforms_to_peak_correctness():
    # delta functions at different positions across waveforms
    n_wfs, n_samples = 5, 200
    peak_positions = [20, 50, 100, 150, 180]
    wf = np.zeros((n_wfs, n_samples))
    for i, p in enumerate(peak_positions):
        wf[i, p] = 1.0

    alignment_idx = 100

    shifted, returned_peaks = psl.align_waveforms_to_peak(
        wf, alignment_idx=alignment_idx, nsamples_output_current_wfs=n_samples
    )

    # returned peak indices must match the input positions
    np.testing.assert_array_equal(returned_peaks, peak_positions)

    # for waveforms whose shift doesn't push the peak outside the output window,
    # the peak must land exactly at alignment_idx after shifting
    for i, p in enumerate(peak_positions):
        shift = alignment_idx - p
        dst = p + shift  # == alignment_idx
        if 0 <= dst < n_samples:
            assert np.argmax(shifted[i]) == alignment_idx, (
                f"waveform {i}: peak expected at {alignment_idx}, "
                f"got {np.argmax(shifted[i])}"
            )


def test_align_waveforms_to_peak_clipping():
    # peak near the start: left-clipping
    wf_left = np.zeros((1, 100))
    wf_left[0, 5] = 1.0
    shifted_left, _ = psl.align_waveforms_to_peak(
        wf_left, alignment_idx=50, nsamples_output_current_wfs=100
    )
    # peak lands at 50; samples originally before index 0 are zero-padded
    assert shifted_left[0, 50] == 1.0
    assert np.sum(shifted_left[0, :50]) == 0.0

    # peak near the end: right-clipping
    wf_right = np.zeros((1, 100))
    wf_right[0, 95] = 1.0
    shifted_right, _ = psl.align_waveforms_to_peak(
        wf_right, alignment_idx=10, nsamples_output_current_wfs=100
    )
    # peak lands at 10; samples originally after index 99 are lost (zero-padded)
    assert shifted_right[0, 10] == 1.0
    assert np.sum(shifted_right[0, 11:]) == 0.0


def test_align_waveforms_to_peak_different_output_length():
    # output shorter than input: exercises length-asymmetric index arithmetic
    n_wfs, n_samples_in, n_samples_out = 3, 200, 150
    # peaks at: 10 (large positive shift → left-clip into shorter output),
    #           80 (fits comfortably),
    #           190 (large negative shift → right-clip into shorter output)
    peak_positions = [10, 80, 190]
    alignment_idx = 75

    wf = np.zeros((n_wfs, n_samples_in))
    for i, p in enumerate(peak_positions):
        wf[i, p] = 1.0

    shifted, returned_peaks = psl.align_waveforms_to_peak(
        wf, alignment_idx=alignment_idx, nsamples_output_current_wfs=n_samples_out
    )

    assert shifted.shape == (n_wfs, n_samples_out)
    np.testing.assert_array_equal(returned_peaks, peak_positions)

    # waveform 1 (peak at 80): shift = 75 - 80 = -5, fits entirely in output
    assert np.argmax(shifted[1]) == alignment_idx
    assert shifted[1, alignment_idx] == 1.0

    # waveform 0 (peak at 10): shift = +65, peak lands at 75
    assert shifted[0, alignment_idx] == 1.0
    assert np.sum(shifted[0, :alignment_idx]) == 0.0  # left side zero-padded

    # waveform 2 (peak at 190): shift = 75 - 190 = -115, dst = 75
    # but src start = 115, src end = min(200, 150+115) = 200 → only 85 samples copied
    assert shifted[2, alignment_idx] == 1.0
    assert np.sum(shifted[2, alignment_idx + 1 :]) == 0.0  # right side zero-padded


def test_align_waveforms_to_peak_awkward_input():
    wf_np = np.zeros((4, 100))
    wf_np[:, 40] = 1.0
    wf_ak = ak.Array(wf_np.tolist())

    shifted_np, peaks_np = psl.align_waveforms_to_peak(
        wf_np, alignment_idx=50, nsamples_output_current_wfs=100
    )
    shifted_ak, peaks_ak = psl.align_waveforms_to_peak(
        wf_ak, alignment_idx=50, nsamples_output_current_wfs=100
    )

    np.testing.assert_array_equal(peaks_np, peaks_ak)
    np.testing.assert_array_almost_equal(shifted_np, shifted_ak)


def test_make_realistic_pulse_shape_lib():
    rng = np.random.default_rng()
    ideal_psl = {
        "r": Array([0.1, 0.2, 0.3], attrs={"units": "m"}),
        "z": Array([-0.1, 0, 0.1], attrs={"units": "m"}),
        "waveform_0": Array(rng.random((3, 1000))),
        "dt": Scalar(10, attrs={"units": "ns"}),
    }

    kernel = psl.build_electronics_response_kernel(1, 0, 100, 100)

    output = psl.make_realistic_pulse_shape_lib(
        ideal_psl,
        kernel,
        500,
        1000,
        mw_pars={
            "length": 48,
            "num_mw": 3,
            "mw_type": 0,
        },
    )
    assert isinstance(output, dict)

    for item in output.values():
        assert isinstance(item, (Array, Scalar))

    assert "r" in output
    assert "z" in output

    assert "waveform_0" in output

    assert output["r"].view_as("np").shape == (3,)
    assert output["z"].view_as("np").shape == (3,)

    assert output["waveform_0"].view_as("np").shape == (3, 1000)
    assert "drift_time_0" in output
    assert output["drift_time_0"].view_as("np").shape == (3,)


_MW_PARS = {"length": 48, "num_mw": 3, "mw_type": 0}


def test_make_realistic_pulse_shape_lib_nan_propagation():
    # Build a 2D waveform map (n_r=3, n_z=4) where one pixel is NaN
    rng = np.random.default_rng(0)
    wfs = rng.random((3, 4, 1000))
    nan_pixel = (1, 2)
    wfs[nan_pixel] = np.nan

    ideal_psl = {
        "r": Array([0.1, 0.2, 0.3], attrs={"units": "m"}),
        "z": Array([-0.15, -0.05, 0.05, 0.15], attrs={"units": "m"}),
        "waveform_0": Array(wfs),
        "dt": Scalar(10, attrs={"units": "ns"}),
    }

    kernel = psl.build_electronics_response_kernel(1, 0, 100, 100)
    output = psl.make_realistic_pulse_shape_lib(
        ideal_psl, kernel, 500, 1000, mw_pars=_MW_PARS
    )

    out_wfs = output["waveform_0"].view_as("np")
    out_dt = output["drift_time_0"].view_as("np")

    # NaN pixel must remain NaN in both outputs
    assert np.all(np.isnan(out_wfs[nan_pixel])), "NaN pixel not propagated to waveform"
    assert np.isnan(out_dt[nan_pixel]), "NaN pixel not propagated to drift time"

    # all other pixels must be finite
    mask = np.zeros((3, 4), dtype=bool)
    mask[nan_pixel] = True
    assert np.all(np.isfinite(out_wfs[~mask])), (
        "unexpected NaN in valid waveform pixels"
    )
    assert np.all(np.isfinite(out_dt[~mask])), (
        "unexpected NaN in valid drift-time pixels"
    )


def test_make_realistic_pulse_shape_lib_3d():
    # Simulate the real input shape coming from Julia via HDF5:
    # Julia writes (n_samples, n_z, n_r); Python reads back (n_r, n_z, n_samples)
    n_r, n_z, n_samples = 4, 5, 1000
    rng = np.random.default_rng(1)
    wfs = rng.random((n_r, n_z, n_samples))

    ideal_psl = {
        "r": Array(np.linspace(0.0, 0.03, n_r), attrs={"units": "m"}),
        "z": Array(np.linspace(0.0, 0.04, n_z), attrs={"units": "m"}),
        "waveform_0": Array(wfs),
        "dt": Scalar(10, attrs={"units": "ns"}),
    }

    kernel = psl.build_electronics_response_kernel(1, 0, 100, 100)
    output = psl.make_realistic_pulse_shape_lib(
        ideal_psl, kernel, 500, 1000, mw_pars=_MW_PARS
    )

    assert output["waveform_0"].view_as("np").shape == (n_r, n_z, n_samples)
    assert output["drift_time_0"].view_as("np").shape == (n_r, n_z)
