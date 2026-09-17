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
import warnings
from collections.abc import Mapping
from pathlib import Path

import awkward as ak
import dbetto
import hist
import lh5
import matplotlib.pyplot as plt
import numpy as np
from dspeed.processors import moving_window_multi
from lgdo import Array, Scalar
from matplotlib.figure import Figure
from numpy.typing import DTypeLike
from reboost import units
from scipy.signal import convolve, fftconvolve

log = logging.getLogger(__name__)


def validate_ssd_scan_grid(file: str, detector: str) -> bool:
    """Validate the structure of an SSD scan file.

    These files can contain PSLs, electronics model parameters or other
    observables and have a common structure::

        /
        └── DETECTOR · struct{grid_info,psl_scan}
            ├── grid_info · struct{dep_min,dep_step,slope_min,slope_step}
            │   ├── dep_min · real
            │   ├── dep_step · real
            │   ├── slope_min · real
            │   └── slope_step · real
            └── psl_scan · struct{slope_0, slope_1, ..., slope_M}
                ├── slope_0 · struct{dep_0, dep_1, ..., dep_N}
                │   ├── dep_0 · struct{...}
                │   ├── dep_1 · struct{...}
                │   :
                │   └── dep_N · struct{...}
                :
                └── slope_M · struct{dep_0, dep_1, ..., dep_N}

    The top-level group is named after the detector. A file can hold more than
    one detector, each validated on its own.

    This structure can either be implemented in LH5 files or in the text
    formats read by :func:`dbetto.utils.load_dict`, i.e. YAML and JSON. In the
    case of LH5 files the structure is implemented as above, in the case of
    text files as a nested dictionary.

    This represents a 2D scan of the simulation over `M + 1` slopes and `N + 1`
    depletion voltage parameters, with an arbitrary data object for each
    combination. This format only defines the grid structure: the only
    requirement on the underlying data structs is that they should all have the
    same structure.

    The groups are numbered from 0 in scan order, so the values behind `slope_i`
    and `dep_j` are `slope_min + i * slope_step` and `dep_min + j * dep_step`,
    with `grid_info` giving the grid in physical units:

    - `dep_min` and `dep_step` are depletion voltages in V. They are stored as
      absolute values, while the scan settings that generate them are given as
      shifts relative to the operational voltage of the detector, see
      :ref:`psl-scan-settings-meta`.
    - `slope_min` and `slope_step` are dimensionless scaling factors applied to
      the non-constant part of the impurity profile.

    Parameters
    ----------
    file
        Path to the LH5, YAML or JSON file to validate.
    detector
        Name of the detector to validate in the file.
    """
    suffix = Path(file).suffix
    # the text formats dbetto can read, i.e. YAML and JSON
    text_suffixes = {
        ext for exts in dbetto.utils.__file_extensions__.values() for ext in exts
    }

    if suffix == ".lh5":

        def list_fields(group: str) -> list[str]:
            return lh5.ls(file, f"{group}/")

    elif suffix in text_suffixes:
        contents = dbetto.utils.load_dict(file)

        def list_fields(group: str) -> list[str]:
            node = contents
            for part in group.split("/"):
                if not isinstance(node, Mapping):
                    return []
                node = node.get(part, {})
            return (
                [f"{group}/{key}" for key in node] if isinstance(node, Mapping) else []
            )

    else:
        msg = (
            f"Cannot validate '{file}': expected a .lh5 file or one of "
            f"{sorted(text_suffixes)}"
        )
        raise ValueError(msg)

    fields = list_fields(detector)

    if f"{detector}/grid_info" not in fields:
        msg = f"Missing 'grid_info' group in {detector} of {file}"
        log.info(msg)
        return False

    grid_info_fields = {f.split("/")[-1] for f in list_fields(f"{detector}/grid_info")}
    required_grid_info_fields = {"dep_min", "dep_step", "slope_min", "slope_step"}
    missing_fields = required_grid_info_fields - grid_info_fields

    if missing_fields:
        msg = f"Missing fields in 'grid_info' of {detector} in {file}: {missing_fields}"
        log.info(msg)
        return False

    if f"{detector}/psl_scan" not in fields:
        msg = f"Missing 'psl_scan' group in {detector} of {file}"
        log.info(msg)
        return False

    slope_fields = list_fields(f"{detector}/psl_scan")
    if not slope_fields:
        msg = f"No slope groups in 'psl_scan' of {detector} in {file}"
        log.info(msg)
        return False

    dep_fields = [
        [d.split("/")[-1] for d in list_fields(slope)] for slope in slope_fields
    ]

    if not all(dep and dep == dep_fields[0] for dep in dep_fields):
        msg = f"Missing 'dep' groups in some slopes of {detector} in {file}"
        log.info(msg)
        return False

    return True


def compare_psl_scans(file1: str, file2: str, detector: str) -> bool:
    """Compare two pulse-shape library scan files.

    This function compares the structure and content of two LH5 files
    containing pulse-shape library scans. It checks for the presence of the
    same detector, the same slope and depletion voltage parameters.

    Parameters
    ----------
    file1
        Path to the first LH5 file to compare.
    file2
        Path to the second LH5 file to compare.
    detector
        Name of the detector to compare in both LH5 files.
    """
    raise NotImplementedError



def load_ideal_psl_scan(psl_file:str)->tuple[dict[str, dict[str, np.ndarray]], dict[str, np.ndarray]]:
    """Load the ideal pulse-shape library scan. """

    dets = lh5.ls(psl_file, "/")

    assert len(dets) == 1
    det = dets[0]

    output = {}

    slopes = lh5.ls(psl_file, f"{det}/psl_scan/")

    for slope_group in slopes:

        slope = slope_group.split("/")[-1]

        depv_groups = lh5.ls(psl_file, f"{slope_group}/")
        output[slope] = {}
        
        for depv in depv_groups:
            depv_name = depv.split("/")[-1]
            output[slope][depv_name] = lh5.read(psl_file, f"{det}/psl_scan/{slope}")

    info = lh5.read(f"{det}/info", psl_file)

    return output, info

def convolve_elecmod_scan(ideal_psls:dict, sigma:float,tau:float,alignment_idx = 1000.,n_samples = 4001,mw_pars = {"length": 48, "num_mw": 3, "mw_type": 0}, dt_data =16):
    """ Convolve ideal PSLs from scan, with the electronics model."""

    output = {}
    for slope, depv_psls in ideal_psls.items():
        psls[slope] = {}
        
        for depv, ideal_psl in depv_psls.items():

            dt = ideal_psl["dt"].value * units.units_convfact(ideal_psl["dt"], "ns")

            rf_kernel = psl.build_electronics_response_kernel(
                dt,
                mu_bandwidth=0,
                sigma_bandwidth=elecmod_,
                tau_rc=tau,
                kernel_start=-100,
            )

            realistic_dict = make_realistic_pulse_shape_lib(
                ideal_map_obj,
                rf_kernel,
                alignment_idx,
                n_samples,
                mw_pars=mw_pars,
                dt_data=dt_data,
                dtype=np.float32,
                kernel_t0_idx=-2 * kernel_start,
            )
            h_aoe, mean_aoe = get_avg_aoe(
                [realistic_dict[k] for k in realistic_dict if "waveform" in k]
            )

            for key in realistic_dict:
                if "waveform" in key:
                    realistic_dict[key] = realistic_dict[key].view_as("np") / mean_aoe
                    
            output[slope][depv] = realistic_dict

   return output


def get_avg_aoe(waveforms: list[np.ndarray]) -> tuple[hist.Hist, float]:
    """Estimate the average A/E from the PSL.

    Estimated as the mode of the distribution of
    the maximum amplitude of each waveform.

    Parameters
    ----------
    waveforms
        List of 3D array of waveforms with shape (n_r, n_z, n_samples)

    Returns
    -------
    hist_aoe
        Histogram of the maximum amplitude distribution.
    avg_aoe
        The average A/E value estimated from the waveforms
    """
    aoe = np.concatenate([np.max(waveform, axis=2).ravel() for waveform in waveforms])
    aoe = aoe[~np.isnan(aoe)]

    amin, amax = float(np.min(aoe)), float(np.max(aoe))
    if amax <= amin:
        # degenerate distribution (e.g. constant amplitudes): pad to a finite
        # range so the histogram axis is valid
        pad = abs(amin) * 1e-3 or 1.0
        amin, amax = amin - pad, amax + pad

    hist_aoe = hist.new.Reg(1000, amin, amax, name="A/E", label="A/E").Double()

    hist_aoe.fill(aoe)
    counts, bin_edges = hist_aoe.to_numpy()
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    # plain float: an np.float64 would upcast the float32 waveforms it normalizes
    avg_aoe = float(bin_centers[np.argmax(counts)])

    return hist_aoe, avg_aoe


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

    ``t=0`` sits at index ``-2 * kernel_start`` (``-kernel_start`` if
    `gaussian_only`), not at the kernel maximum (the Gaussian makes the kernel
    non-causal).

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
    nsamples_output_wfs: int,
    *,
    peak_indices: np.ndarray | None = None,
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
    nsamples_output_wfs
        The total length of the resulting aligned current waveforms
    peak_indices
        If not `None` use this as the indices for alignment.

    Returns
    -------
    shifted_wfs
        2D array of shifted waveforms
    peak_indices
        1D array containing the original peak index for each current waveform

    """
    # Convert to NumPy for matrix operations
    wfs = ak.to_numpy(wf_input) if isinstance(wf_input, ak.Array) else wf_input

    n_wfs, n_samples = wfs.shape

    if peak_indices is None:
        peak_indices = np.argmax(wfs, axis=1)

    shifts = alignment_idx - peak_indices  # (n_wfs,)

    # For each destination index d, the corresponding source index is d - shift[i].
    # Build a (n_wfs, n_samples_out) source-index matrix, mask out-of-bounds positions
    # with 0 (safe to index), then zero out those positions in the result.
    d = np.arange(nsamples_output_wfs)  # (n_samples_out,)
    src_idx = d[np.newaxis, :] - shifts[:, np.newaxis]  # (n_wfs, n_samples_out)
    valid = (src_idx >= 0) & (src_idx < n_samples)
    src_idx_safe = np.clip(src_idx, 0, n_samples - 1)

    shifted_wfs = np.where(
        valid, wfs[np.arange(n_wfs)[:, np.newaxis], src_idx_safe], 0.0
    )

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


def process_ideal_waveforms(
    wfs: np.ndarray,
    rf_kernel: np.ndarray,
    dt: float,
    alignment_idx: int,
    nsamples_output: int,
    mw_pars: Mapping[str, int] | None = None,
    dt_data: float = 16.0,
    return_mode: str = "current",
    dtype: DTypeLike = np.float32,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply electronics response and DSP chain to ideal charge waveforms.

    Convolve -> differentiate -> MWA -> align to peak.

    Parameters
    ----------
    wfs
        Charge waveforms, shape ``(n_wfs, n_samples)``.
    rf_kernel
        System response kernel from
        :func:`build_electronics_response_kernel`.
    dt
        Time step of the ideal waveforms in ns.
    alignment_idx
        Sample index where the current peak is placed after alignment.
    nsamples_output
        Length of the output waveforms.
    mw_pars
        MWA parameters: ``length``, ``num_mw``, ``mw_type``. Defaults to the
        settings of the LEGEND-200 production DSP chain.
    dt_data
        Data sampling time step in ns. Defaults to the LEGEND-200 digitizer
        period.
    return_mode
        Whether to extract the "current" or the "charge" waveform.
    dtype
        Floating-point type used throughout the processing chain. float32 is
        more than enough for the ~1% A/E it feeds and halves the library
        footprint.

    Returns
    -------
    aligned_currents
        Aligned current waveforms, shape ``(n_wfs, nsamples_output)``.
    current_peak_indices
        Index of the current peak for each waveform after MWA and before
        alignment, shape ``(n_wfs,)``.

    """
    if mw_pars is None:
        mw_pars = {"length": 48, "num_mw": 3, "mw_type": 0}

    convolved = ak.to_numpy(
        apply_electronics_response(
            np.asarray(wfs, dtype=dtype), rf_kernel.astype(dtype, copy=False)
        )
    )

    # Derivative (charge -> current), scaled to data sampling units
    current = np.diff(convolved, axis=-1, prepend=0) * (dt_data / dt)

    # Moving window average
    mwa_out = np.zeros_like(current, dtype=dtype)
    moving_window_multi(
        current.astype(dtype, copy=False),
        mw_pars["length"],
        mw_pars["num_mw"],
        mw_pars["mw_type"],
        mwa_out,
    )

    # Record peak indices after MWA: this is both the reference used for the
    # alignment below and the analogue of tp_aoe_max, which the production DSP
    # chain extracts from the moving-window-averaged current (curr_av)
    current_peak_indices = np.argmax(mwa_out, axis=1)

    aligned_currents, _ = align_waveforms_to_peak(
        mwa_out, alignment_idx, nsamples_output, peak_indices=current_peak_indices
    )
    aligned_charges, _ = align_waveforms_to_peak(
        convolved,
        alignment_idx,
        nsamples_output,
        peak_indices=current_peak_indices,
    )
    if return_mode == "current":
        return aligned_currents, current_peak_indices
    if return_mode == "charge":
        return aligned_charges, current_peak_indices
    msg = f"Invalid return_mode '{return_mode}', expected 'current' or 'charge'"
    raise ValueError(msg)


def make_realistic_pulse_shape_lib(
    ideal_pulse_shape_lib_obj: Mapping[str, Array | Scalar],
    rf_kernel: np.ndarray,
    alignment_idx: int,
    nsamples_output_current_wfs: int,
    mw_pars: Mapping[str, int] | None = None,
    dt_data: float = 16.0,
    dtype: DTypeLike = np.float32,
    *,
    kernel_t0_idx: int,
) -> dict[str, Array | Scalar]:
    """Apply the waveform post-processing chain to generate a realistic waveform map.

    Starts from an ideal waveform map and performs the following steps:

    1. Converts coordinates (m to mm)
    2. Convolves with system response
    3. Aligns by Peak Time
    4. Calculates the drift time: filtered-current peak (the alignment sample)
       measured from the energy deposition time, i.e. minus `kernel_t0_idx`

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
        Parameters for the moving window average step, with keys:

        - length: The length of the moving window (in samples)
        - num_mw: The number of moving windows to use in the moving window
          average
        - mw_type: The type of moving window to apply (see
          ``dspeed.processors.moving_window_multi`` for details)

        Defaults to the settings of the LEGEND-200 production DSP chain.
    dt_data
        The time step of the original data waveforms (in ns), used to scale
        the derivative. Defaults to the LEGEND-200 digitizer period.
    dtype
        Floating-point type of the waveform and drift-time samples, both in
        memory and in the output library. float32 is more than enough for the
        ~1% A/E it feeds and halves the library footprint.
    kernel_t0_idx
        Index of ``t=0`` in `rf_kernel`, ``-2 * kernel_start`` for
        :func:`build_electronics_response_kernel`.

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

    for coord in ["r", "z"]:
        if coord in ideal_pulse_shape_lib_obj:
            realistic_pulse_shape_lib[coord] = Array(
                units.units_conv_ak(ideal_pulse_shape_lib_obj[coord], "mm").to_numpy(),
                attrs={"units": "mm"},
            )

    # Search for waveform keys (general - I can have any angle in the ideal map)
    keys_to_convolve = [k for k in ideal_pulse_shape_lib_obj if "waveform" in k]

    for key in keys_to_convolve:
        log.info("Processing %s...", key)

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

        # Full processing chain: convolve -> differentiate -> MWA -> align
        curr_aligned, current_peak_indices = process_ideal_waveforms(
            wfs_flat,
            rf_kernel,
            dt,
            alignment_idx,
            nsamples_output_current_wfs,
            mw_pars,
            dt_data,
            dtype=dtype,
        )

        drift_indices = current_peak_indices - kernel_t0_idx
        drift_times_flat = drift_indices * dt
        drift_times_2d = drift_times_flat.reshape(original_shape[:-1]).astype(dtype)

        # Restore NaN for invalid pixels in drift time
        drift_times_2d[nan_mask] = np.nan
        dt_key = key.replace("waveform", "drift_time")
        realistic_pulse_shape_lib[dt_key] = Array(drift_times_2d, attrs={"units": "ns"})

        new_length = curr_aligned.shape[-1]
        new_shape = (*original_shape[:-1], new_length)
        wfs_out = curr_aligned.reshape(new_shape).astype(dtype, copy=False)

        # Restore NaN for invalid pixels in waveforms
        wfs_out[nan_mask] = np.nan
        realistic_pulse_shape_lib[key] = Array(wfs_out, attrs={"units": ""})

    _check_pulse_shape_lib_keys(realistic_pulse_shape_lib)
    return realistic_pulse_shape_lib


# Plots


def plot_rz_scan(
    pulse_shape_lib: Mapping[str, Array | Scalar],
    angle_deg: int,
    detector_id: str,
    scan: str = "r",
    step: int = 1,
    xlim: tuple[float, float] | None = None,
) -> tuple[Figure, plt.Axes]:
    """Plot an R or Z scan of waveforms from a pulse-shape library.

    For a given azimuthal angle, produces one figure scanning over
    radial or axial positions at a fixed index on the other axis.
    Waveforms are color-coded by spatial coordinate.

    Parameters
    ----------
    pulse_shape_lib
        Mapping containing the pulse-shape library (ideal or realistic),
        as returned by :func:`make_realistic_pulse_shape_lib` or read from
        the ideal LH5 file.  Must contain ``r``, ``z``, ``dt`` keys and at
        least one ``waveform_<angle>_deg`` key.
    angle_deg
        Azimuthal angle in degrees, used to select the
        ``waveform_<angle>_deg`` key.
    detector_id
        Detector name, e.g. ``"V03422A"``.
    scan
        ``"r"`` to scan over radial positions at fixed Z,
        ``"z"`` to scan over axial positions at fixed R.
    step
        Stride for subsampling spatial positions (default: plot every position).
    xlim
        x-axis limits in ns. If ``None``, matplotlib auto-scales.


    Returns
    -------
    fig : Figure
    ax : Axes
    """
    import matplotlib.colors as mcolors  # noqa: PLC0415
    from matplotlib import cm  # noqa: PLC0415

    if scan not in ("r", "z"):
        msg = f"scan must be 'r' or 'z', got {scan!r}"
        raise ValueError(msg)

    wf_key = f"waveform_{str(angle_deg).zfill(3)}_deg"
    if wf_key not in pulse_shape_lib:
        msg = f"Key '{wf_key}' not found in pulse-shape library"
        raise KeyError(msg)

    wfs = pulse_shape_lib[wf_key].nda
    r_mm = pulse_shape_lib["r"].nda
    z_mm = pulse_shape_lib["z"].nda
    dt = pulse_shape_lib["dt"].value
    t0 = pulse_shape_lib["t0"].value if "t0" in pulse_shape_lib else 0.0

    Nr, Nz, Nt = wfs.shape
    time_axis = t0 + np.arange(Nt) * dt

    if scan == "r":
        idx_fix = Nz // 2
        axis_vals = r_mm
        wf_slice = wfs[:, idx_fix, :]
        label = "Radius [mm]"
        title_ctx = f"R scan - Fixed Z = {z_mm[idx_fix]:.1f} mm"
    else:
        idx_fix = Nr // 4
        axis_vals = z_mm
        wf_slice = wfs[idx_fix, :, :]
        label = "Z [mm]"
        title_ctx = f"Z scan - Fixed R = {r_mm[idx_fix]:.1f} mm"

    cmap = "magma_r"
    norm = mcolors.Normalize(vmin=np.min(axis_vals), vmax=np.max(axis_vals))
    cmap_obj = plt.get_cmap(cmap)

    fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")

    for i in range(0, len(axis_vals), step):
        color = cmap_obj(norm(axis_vals[i]))
        ax.plot(time_axis, wf_slice[i], color=color, lw=1.5, linestyle="-")

    if xlim is not None:
        ax.set_xlim(xlim)
    ax.set_title(f"{detector_id} - {title_ctx} ({angle_deg}°)")
    ax.set_xlabel("Time [ns]")
    ax.set_ylabel("Signal [A.U.]")
    ax.grid(visible=True, which="both", axis="both", linestyle="--", alpha=0.5)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(label)

    return fig, ax


def symmetrize(a: np.ndarray) -> np.ndarray:
    """Mirror a half ``[n_r, n_z]`` R/Z map about ``r = 0`` into a full image."""
    a = a.T
    return np.concatenate((np.fliplr(a), a), axis=1)


def plot_aoe_rz_map(
    pulse_shape_lib: Mapping[str, Array | Scalar],
    detector_id: str,
    *,
    hpge_profile: object | None = None,
) -> Figure:
    """Plot the A/E R/Z map(s) of a pulse-shape library.

    For each azimuthal angle present in the library, computes the per-pixel A/E
    as the maximum of the (energy-normalized) current waveform over the time
    axis and renders it as a symmetrized R/Z heatmap, styled like the HPGe
    drift-time-map validation plot. When both the ``<100>`` (0 deg) and
    ``<110>`` (45 deg) angles are present, an additional ``<100>/<110>`` ratio
    panel is appended.

    Parameters
    ----------
    pulse_shape_lib
        Mapping containing the pulse-shape library, as returned by
        :func:`make_realistic_pulse_shape_lib`. Must contain ``r``, ``z`` (in
        mm) and at least one ``waveform_<angle>_deg`` key.
    detector_id
        Detector name, e.g. ``"V03422A"``, used in the panel titles.
    hpge_profile
        Optional detector geometry object (as returned by
        ``pygeomhpges.make_hpge``). When given, its profile is overlaid on every
        panel.

    Returns
    -------
    fig : Figure
    """
    import pygeomhpges.draw  # noqa: PLC0415

    axis_labels = {0: "<100>", 45: "<110>"}

    angles = sorted(int(k.split("_")[1]) for k in pulse_shape_lib if "waveform" in k)
    if not angles:
        msg = "no 'waveform_<angle>_deg' keys found in pulse-shape library"
        raise KeyError(msg)

    r = pulse_shape_lib["r"].nda  # already in mm
    z = pulse_shape_lib["z"].nda

    # per-angle A/E maps (max current amplitude over the time axis)
    images = {}
    for angle in angles:
        wfs = pulse_shape_lib[f"waveform_{angle:03d}_deg"].nda
        with warnings.catch_warnings():
            # invalid pixels are all-NaN waveforms; keep them NaN (rendered blank)
            warnings.simplefilter("ignore", category=RuntimeWarning)
            aoe = np.nanmax(wfs, axis=-1)
        images[angle] = symmetrize(aoe)

    has_ratio = {0, 45}.issubset(angles)
    n_panels = len(angles) + (1 if has_ratio else 0)

    extent = (-r.max(), r.max(), z.min(), z.max())
    stacked = list(images.values())
    vmin, vmax = np.nanmin(stacked), np.nanmax(stacked)

    fig, axes = plt.subplots(
        ncols=n_panels,
        figsize=(4 * n_panels, 4),
        sharey=True,
        squeeze=False,
    )
    axes = axes[0]

    def plot(ax, img, title, *, cmap, vmin=None, vmax=None):
        im = ax.imshow(
            img,
            origin="lower",
            extent=extent,
            aspect="equal",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        if hpge_profile is not None:
            pygeomhpges.draw.plot_profile(
                hpge_profile,
                axes=ax,
                marker=None,
                linewidth=1,
                color="black",
            )

        xmin, xmax, ymin, ymax = extent

        f = 0.04
        ax.set_xlim(xmin - f * (xmax - xmin), xmax + f * (xmax - xmin))
        ax.set_ylim(ymin - f * (ymax - ymin), ymax + f * (ymax - ymin))

        ax.set_xlabel("r (mm)")
        ax.set_title(f"{detector_id} · {title}")

        return im

    im = None
    for ax, angle in zip(axes[: len(angles)], angles, strict=True):
        im = plot(
            ax,
            images[angle],
            axis_labels.get(angle, f"{angle}°"),
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
        )

    axes[0].set_ylabel("z (mm)")

    # shared A/E colorbar across all angle panels
    fig.colorbar(im, ax=axes[: len(angles)], label="A/E")

    if has_ratio:
        ratio = np.divide(
            images[0],
            images[45],
            out=np.full_like(images[0], np.nan),
            where=images[45] > 0,
        )
        im_ratio = plot(axes[-1], ratio, "<100> / <110>", cmap="coolwarm")
        fig.colorbar(im_ratio, ax=axes[-1], label="ratio")

    return fig
