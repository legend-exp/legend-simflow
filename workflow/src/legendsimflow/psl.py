# Copyright (C) 2026 Luigi Pertoldi <gipert@pm.me>
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
from collections.abc import Mapping

import awkward as ak
import numpy as np
from lgdo import Array, Scalar
from reboost import units
from scipy.signal import convolve, fftconvolve

logger = logging.getLogger(__name__)


def build_electronics_response_kernel(
    dt: float,
    mu_bandwidth: float,
    sigma_bandwidth: float,
    tau_rc: float,
    gaussian_only: bool = False,
    *,
    kernel_length: int = 600,
    kernel_start: int = -100,
) -> np.ndarray:
    """Create the system response kernel (gaussian + exponential decay).

    This is obtained by convolving a Gaussian (representing the digitizer
    bandwidth) with a causal exponential decay (representing the preamplifier
    response). The kernel is normalized to have a sum of 1.

    Note
    ----
    The 'full' mode of convolution results in a length of 2*kernel_length - 1.
    If `gaussian_only` is True, the kernel will have a length of
    `kernel_length` and will only contain the Gaussian component, since no
    convolution is performed.

    Parameters
    ----------
    dt
        The time step between samples in the waveform (in ns)
    mu_bandwidth
        The mean of the Gaussian representing the digitizer bandwidth (in ns)
    sigma_bandwidth
        The standard deviation of the Gaussian representing the digitizer
        bandwidth (in ns)
    tau_rc
        The time constant of the exponential decay representing the
        preamplifier response (in ns)
    gaussian_only
        If True, only use the Gaussian component (default is False)
    kernel_length
        The total length of the response kernel in samples (default is 600)
    kernel_start
        The starting index of the kernel relative to the waveform (default is
        -100, meaning the kernel will cover from -100 to 500 samples)

    Returns
    -------
    rf
        The normalized response kernel
    """

    # Validate inputs
    if dt <= 0:
        msg = f"dt must be positive, got {dt}"
        raise ValueError(msg)
    if sigma_bandwidth <= 0:
        msg = f"sigma_bandwidth must be positive, got {sigma_bandwidth}"
        raise ValueError(msg)
    if not gaussian_only and tau_rc <= 0:
        msg = f"tau_rc must be positive, got {tau_rc}"
        raise ValueError(msg)

    # Convert to samples
    mu_samples = mu_bandwidth / dt
    sigma_samples = sigma_bandwidth / dt
    tau_samples = tau_rc / dt

    # Define the sample range
    x = np.arange(kernel_start, kernel_start + kernel_length)

    # Compute Gaussian response (digitizer bandwidth)
    rf_digi = np.exp(-0.5 * ((x - mu_samples) / sigma_samples) ** 2)

    # Validate before normalization
    digi_sum = np.sum(rf_digi)
    if digi_sum == 0:
        msg = (
            "Digitizer response normalization failed: sum is zero. "
            "This may indicate numerical underflow with the given sigma_bandwidth."
        )
        raise ValueError(msg)

    rf_digi /= digi_sum

    if gaussian_only:
        return rf_digi

    # Compute causal exponential decay (preamplifier response) - response is zero for x < 0
    rf_preamp = np.zeros_like(x, dtype=float)
    rf_preamp[x >= 0] = np.exp(-x[x >= 0] / tau_samples)

    # Validate before normalization
    preamp_sum = np.sum(rf_preamp)
    if preamp_sum == 0:
        msg = (
            "Preamp response normalization failed: sum is zero. "
            "This may indicate numerical underflow with the given tau_rc."
        )
        raise ValueError(msg)
    rf_preamp /= preamp_sum

    # Convolve preamp and digitizer responses - 'full' mode results in a length of (len(rf_preamp) + len(rf_digi) - 1)
    return convolve(rf_preamp, rf_digi, mode="full")


def apply_electronics_response(
    wf_array: ak.Array | np.ndarray, rf_kernel: np.ndarray, batch_size: int = 50000
) -> ak.Array:
    """Vectorized convolution using FFT with batching to save memory.

    Parameters
    ----------
    wf_array
        Array of waveforms (all of same length)
    rf_kernel
        The response kernel (gaussian + exponential)
    batch_size
        Number of waveforms to process at once (default is 50,000)

    Returns
    -------
    convolved_wfs
        The convolved waveforms as an Awkward Array
    """
    n_events = len(wf_array)
    # Reshape kernel for broadcasting (1, Kernel_Length)
    kernel_2d = rf_kernel[np.newaxis, :]

    convolved_results = []

    # Process in batches to prevent memory overflow
    for i in range(0, n_events, batch_size):
        batch = wf_array[i : i + batch_size]
        wf_matrix = np.asarray(batch)
        wf_conv = fftconvolve(wf_matrix, kernel_2d, mode="full")

        # Slice to original length to maintain timing relative to start. This preserves the initial points to emulate baseline
        wf_conv_sliced = wf_conv[:, : wf_matrix.shape[1]]
        convolved_results.append(wf_conv_sliced)

    # Concatenate batches and return as Awkward Array
    return ak.concatenate(convolved_results)


def align_waveforms_to_peak(
    wf_input: ak.Array | np.ndarray,
    alignment_idx: int,
    nsamples_output_current_wfs: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Align an array of waveforms by shifting their maximum to a fixed index.

    No normalization is performed; raw amplitudes are preserved.

    Note
    ----
    The output peak_indices is not the original drift time as the current
    waveform inherits baseline from convolution.

    Parameters
    ----------
    wf_input
        Input array of current waveforms
    alignment_idx
        The index in the output array where the peak will be placed
    nsamples_output_current_wfs
        The total length of the resulting aligned current waveforms

    Returns
    -------
    shifted_wfs
        2D array of shifted waveforms
    peak_indices
        1D array containing the original peak index for each current waveform
    """

    # Convert to NumPy for matrix operations
    wfs = ak.to_numpy(wf_input) if isinstance(wf_input, ak.Array) else wf_input

    # Output Array - Indices not covered by the shifted signal are automatically padded with 0.0
    n_wfs, n_samples = wfs.shape
    shifted_wfs = np.zeros((n_wfs, nsamples_output_current_wfs), dtype=float)

    # Find Amax
    peak_indices = np.argmax(wfs, axis=1)

    # Shift each waveform
    for i in range(n_wfs):
        peak = peak_indices[i]
        shift = alignment_idx - peak

        # Bounds to prevent IndexErrors
        start_src = max(0, -shift)
        end_src = min(n_samples, nsamples_output_current_wfs - shift)

        if start_src < end_src:
            start_dst = start_src + shift
            end_dst = end_src + shift

            shifted_wfs[i, start_dst:end_dst] = wfs[i, start_src:end_src]

    return shifted_wfs, peak_indices


def _check_pulse_shape_lib_keys(
    pulse_shape_lib: Mapping[str, Array | Scalar],
) -> None:
    """Validate that the waveform map contains the required keys with correct types."""
    required_keys = {"r", "z", "dt"}
    missing_keys = required_keys - set(pulse_shape_lib.keys())

    if missing_keys:
        msg = f"Ideal waveform map is missing required keys: {missing_keys}"
        raise ValueError(msg)

    for key in required_keys:
        if not isinstance(pulse_shape_lib[key], (Array, Scalar)):
            msg = f"Key '{key}' must be of type Array or Scalar, got {type(pulse_shape_lib[key])}"
            raise TypeError(msg)

    for key, item in pulse_shape_lib.items():
        if not isinstance(item, (Array, Scalar)):
            msg = f"Key '{key}' must be of type Array or Scalar, got {type(item)}"
            raise TypeError(msg)


def make_realistic_pulse_shape_lib(
    ideal_pulse_shape_lib_obj: Mapping[str, Array | Scalar],
    rf_kernel: np.ndarray,
    alignment_idx: int,
    nsamples_output_current_wfs: int,
    mw_pars: dict[str, float | int],
    dt_data: float = 1.0,
) -> dict[str, Array | Scalar]:
    """Apply the waveform post-processing chain to generate a realistic waveform map.

    Starts from an ideal waveform map and performs the following steps:

    1. Converts coordinates (m to mm)
    2. Convolves with system response
    3. Aligns by Peak Time
    4. Calculates compensated Drift Time

    Parameters
    ----------
    ideal_pulse_shape_lib_obj
        Mapping containing the ideal waveform map with coordinates and waveforms.

        This should have the following format:

        - r: 1D array of radial coordinates
        - z: 1D array of axial coordinates
        - dt: Time step between samples in the waveforms
        - waveform_X: 3D array of ideal charge waveforms for angle X
          (shape: [n_z, n_r, n_samples])
    rf_kernel
        The system response kernel (from :func:`build_electronics_response_kernel`)
    alignment_idx
        The index in the output array where waveform peaks will be aligned
    nsamples_output_current_wfs
        The total length of the resulting aligned current waveforms
    mw_pars
        Dictionary of parameters for the moving window average step, with keys:

        - length: The length of the moving window (in samples)
        - num_mw: The number of moving windows to use in the moving window
          average
        - mw_type: The type of moving window to apply (see
          :func:`dspeed.processors.moving_window_multi` for details)
    dt_data
        The time step of the original data waveforms (in ns), used to scale
        the derivative.

    Returns
    -------
    realistic_pulse_shape_lib
        Struct containing the processed realistic waveform map with the
        following keys:

        - r: 1D array of radial coordinates
        - z: 1D array of axial coordinates
        - t0: Global time offset applied to align waveforms
        - waveform_X: 3D array of processed current waveforms for angle X
          (shape: [n_r, n_z, nsamples_output_current_wfs], spatial axes
          reversed relative to Julia due to HDF5 column-/row-major conversion)
        - drift_time_X: 2D array of calculated drift times for angle X
          (shape: [n_r, n_z])
    """

    _check_pulse_shape_lib_keys(ideal_pulse_shape_lib_obj)

    realistic_pulse_shape_lib = {}
    dt = ideal_pulse_shape_lib_obj["dt"].value * units.units_convfact(
        ideal_pulse_shape_lib_obj["dt"], "ns"
    )

    realistic_pulse_shape_lib["dt"] = Scalar(dt, attrs={"units": "ns"})
    realistic_pulse_shape_lib["t0"] = Scalar(
        -1.0 * alignment_idx * dt, attrs={"units": "ns"}
    )  # Set the global t0 relative to alignment index

    # Delay of the response kernel
    kernel_delay_idx = np.argmax(rf_kernel)

    for coord in ["r", "z"]:
        if coord in ideal_pulse_shape_lib_obj:
            realistic_pulse_shape_lib[coord] = Array(
                units.units_conv_ak(ideal_pulse_shape_lib_obj[coord], "mm").to_numpy(),
                attrs={"units": "mm"},
            )

    # Search for waveform keys (general - I can have any angle in the ideal map)
    keys_to_convolve = [k for k in ideal_pulse_shape_lib_obj if "waveform" in k]

    for key in keys_to_convolve:
        logger.info("Processing %s...", key)

        # Extract and prepare data
        ideal_wfs = ideal_pulse_shape_lib_obj[key]
        ideal_wfs_arr = ideal_wfs.view_as("np")

        # Julia writes arrays in column-major (Fortran) order; HDF5/Python reads them
        # back in row-major (C) order, reversing the axis order. Julia stores waveforms
        # as (n_samples, n_z, n_r), which Python reads as (n_r, n_z, n_samples). The
        # code assumes samples are on the last axis, which is satisfied after this
        # reversal regardless of the spatial axis ordering.

        # Record NaN mask before zeroing (shape: [n_r, n_z], True where pixel is invalid)
        nan_mask = np.isnan(ideal_wfs_arr).any(axis=-1)
        ideal_wfs_arr = np.nan_to_num(ideal_wfs_arr, nan=0.0)

        original_shape = ideal_wfs_arr.shape
        wfs_flat = ideal_wfs_arr.reshape(-1, original_shape[-1])

        # Convolution electronics
        convolved_ak = apply_electronics_response(wfs_flat, rf_kernel)
        wfs_convolved = ak.to_numpy(convolved_ak)

        # Derivative (Charge -> Current)
        wfs_current = np.diff(wfs_convolved, axis=-1, prepend=0) * (dt_data / dt)

        # Calculate drift time from Amax
        current_peak_indices = np.argmax(wfs_current, axis=1)
        drift_indices = current_peak_indices - kernel_delay_idx
        drift_times_flat = drift_indices * dt
        drift_times_2d = drift_times_flat.reshape(original_shape[:-1]).astype(
            np.float64
        )

        # Restore NaN for invalid pixels in drift time
        drift_times_2d[nan_mask] = np.nan
        dt_key = key.replace("waveform", "drift_time")
        realistic_pulse_shape_lib[dt_key] = Array(drift_times_2d, attrs={"units": "ns"})

        # Moving Window Average
        # NOTE: import is intentionally deferred here to work around a bug in dspeed
        # where importing it at module level causes a pint unit registry conflict.
        from dspeed.processors import moving_window_multi  # noqa: PLC0415

        wf_in = wfs_current.astype(float, copy=False)
        wfs_mwa = np.zeros_like(wf_in)
        moving_window_multi(
            wf_in, mw_pars["length"], mw_pars["num_mw"], mw_pars["mw_type"], wfs_mwa
        )

        # Align current waveforms to Amax
        wfs_aligned, _ = align_waveforms_to_peak(
            wfs_mwa, alignment_idx, nsamples_output_current_wfs
        )

        # Reshape
        new_length = wfs_aligned.shape[-1]
        new_shape = (*original_shape[:-1], new_length)
        wfs_out = wfs_aligned.reshape(new_shape).astype(np.float64)

        # Restore NaN for invalid pixels in waveforms
        wfs_out[nan_mask] = np.nan
        realistic_pulse_shape_lib[key] = Array(wfs_out)

    _check_pulse_shape_lib_keys(realistic_pulse_shape_lib)
    return realistic_pulse_shape_lib
