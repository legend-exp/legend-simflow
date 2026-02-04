from __future__ import annotations

import copy
import logging
from pathlib import Path

import awkward as ak
import numpy as np
from dspeed.processors import moving_window_multi
from lgdo import Array, Scalar, Struct, lh5
from scipy.signal import convolve, fftconvolve

logger = logging.getLogger(__name__)


def create_response_kernel(
    dt: float,
    mu_bandwidth: float,
    sigma_bandwidth: float,
    tau_rc: float,
    gaussian_only: bool = False,
) -> np.ndarray:
    """
    Create the system response kernel (gaussian + exponential decay)

    Parameters
    ----------
    dt : float
        The time step between samples in the waveform (in ns)
    mu_bandwidth : float
        The mean of the Gaussian representing the digitizer bandwidth (in ns)
    sigma_bandwidth : float
        The standard deviation of the Gaussian representing the digitizer bandwidth (in ns)
    tau_rc : float
        The time constant of the exponential decay representing the preamplifier response (in ns)
    gaussian_only : bool, optional
        If True, only use the Gaussian component (default is False)

    Returns
    -------
    rf : np.ndarray
        The normalized response kernel
    """

    # Validate inputs
    if dt <= 0:
        raise ValueError(f"dt must be positive, got {dt}")
    if sigma_bandwidth <= 0:
        raise ValueError(f"sigma_bandwidth must be positive, got {sigma_bandwidth}")
    if not gaussian_only and tau_rc <= 0:
        raise ValueError(f"tau_rc must be positive, got {tau_rc}")

    # Convert to samples
    mu_samples = mu_bandwidth / dt
    sigma_samples = sigma_bandwidth / dt
    tau_samples = tau_rc / dt

    # Define the sample range
    x = np.arange(-100, 500)

    # Compute Gaussian response (digitizer bandwidth)
    rf_digi = np.exp(-0.5 * ((x - mu_samples) / sigma_samples) ** 2)
    
    # Validate before normalization
    digi_sum = np.sum(rf_digi)
    if digi_sum == 0:
        raise ValueError(
            "Digitizer response normalization failed: sum is zero. "
            "This may indicate numerical underflow with the given sigma_bandwidth."
        )
    rf_digi /= digi_sum

    if gaussian_only:
        return rf_digi

    # Compute causal exponential decay (preamplifier response) - response is zero for x < 0
    rf_preamp = np.zeros_like(x, dtype=float)
    rf_preamp[x >= 0] = np.exp(-x[x >= 0] / tau_samples)
    
    # Validate before normalization
    preamp_sum = np.sum(rf_preamp)
    if preamp_sum == 0:
        raise ValueError(
            "Preamp response normalization failed: sum is zero. "
            "This may indicate numerical underflow with the given tau_rc."
        )
    rf_preamp /= preamp_sum

    # Convolve preamp and digitizer responses - 'full' mode results in a length of (len(rf_preamp) + len(rf_digi) - 1)
    return convolve(rf_preamp, rf_digi, mode="full")


def convolve_waveforms(
    wf_array: ak.Array, rf_kernel: np.ndarray, batch_size: int = 50000
) -> ak.Array:
    """
    Vectorized convolution using FFT with batching to save memory.

    Parameters
    ----------
    wf_array : ak.Array
        Array of waveforms (all of same length)
    rf_kernel : np.ndarray
        The response kernel (gaussian + exponential)
    batch_size : int, optional
        Number of waveforms to process at once (default is 50,000)

    Returns
    -------
    convolved_wfs : ak.Array
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

        # Slicing
        # Slice to original length to maintain timing relative to start. This preserves the initial points to emulate baseline
        wf_conv_sliced = wf_conv[:, : wf_matrix.shape[1]]
        convolved_results.append(wf_conv_sliced)

    # Concatenate batches and return as Awkward Array
    return ak.concatenate(convolved_results)


def to_current(wf_array: np.ndarray, dt: float, dt_DATA: float = 16) -> np.ndarray:
    """
    Computes the derivative of each charge waveform to convert to current.

    Parameters
    ----------
    wf_array : np.ndarray
        2D array of charge waveforms
    dt : float
        Time step between samples (in ns)
    dt_DATA : float, optional
        Data sampling time for scaling (default is 16)

    Returns
    -------
    derivatives : np.ndarray
        2D array containing the current waveforms - same shape as input
    """

    return (
        np.diff(wf_array, axis=1, prepend=0) / dt
    ) * dt_DATA  # [current] = 1/data_sample


def apply_mwa(wf_array: np.ndarray) -> np.ndarray:
    """
    Apply moving window average to each current waveform in the array.

    Parameters
    ----------
    wf_array : np.ndarray
        2D array of current waveforms [n_events, n_samples]

    Returns
    -------
    mwa_currents : np.ndarray
        Array of filtered waveforms - same shape as input
    """
    wf_in = np.ascontiguousarray(wf_array, dtype=float)

    # Preallocate output array 
    # REVIEW: Does this preallocation cause memory issues with typical dataset sizes?
    # If yes, consider batching approach like in convolve_waveforms()
    mwa_currents = np.zeros_like(wf_in, dtype=float)

    moving_window_multi(
        wf_in, 48, 3, 0, mwa_currents
    ) # REVIEW: Are these MWA parameters appropriate?
    return mwa_currents


def shift_array(
    wf_input: ak.Array | np.ndarray,
    alignment_idx: int,
    nsamples_output_current_wfs: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Aligns an array of waveforms by shifting their maximum to a fixed index.
    No normalization is performed; raw amplitudes are preserved.

    Parameters
    ----------
    wf_input : ak.Array or np.ndarray
        Input array of current waveforms
    alignment_idx : int
        The index in the output array where the peak will be placed
    nsamples_output_current_wfs : int
        The total length of the resulting aligned current waveforms

    Returns
    -------
    shifted_wfs : np.ndarray
        2D array of shifted waveforms
    peak_indices : np.ndarray
        1D array containing the original peak index for each current waveform.
        NOTE: This is NOT the original drift time as current waveform inherits baseline from convolution!
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


def generate_realistic_map(
    ideal_wf_map_obj: Struct,
    rf_kernel: np.ndarray,
    alignment_idx: int,
    nsamples_output_current_wfs: int,
) -> dict:
    """
    Applies the waveform post-processing chain to generate a realistic waveform
    map starting from an ideal one:

    1. Converts coordinates (m to mm) ######################### Correct? Ideal wf map is in meters while reboost wants mm????
    2. Convolves with system response
    3. Aligns by Peak Time
    4. Calculates compensated Drift Time

    Parameters
    ----------
    ideal_wf_map_obj : lgdo.Struct
        LGDO Struct containing the ideal waveform map with coordinates and waveforms
    rf_kernel : np.ndarray
        The system response kernel (from create_response_kernel)
    alignment_idx : int
        The index in the output array where waveform peaks will be aligned
    nsamples_output_current_wfs : int
        The total length of the resulting aligned current waveforms

    Returns
    -------
    realistic_wf_map : dict
        Dictionary containing the processed realistic waveform map with:
        - Converted coordinates (mm)
        - Realistic current waveforms
        - Drift times for each waveform
        - Updated t0 reference
    """

    realistic_wf_map = dict(ideal_wf_map_obj)
    dt = (
        ideal_wf_map_obj["dt"].value
        if hasattr(ideal_wf_map_obj["dt"], "value")
        else ideal_wf_map_obj["dt"]
    )
    realistic_wf_map["t0"] = (
        -1.0 * alignment_idx * dt
    )  # Set the global t0 relative to alignment index

    # Delay of the response kernel
    kernel_delay_idx = np.argmax(rf_kernel)

    # COORDINATE FIX: Convert m to mm ######################### Correct?
    for coord in ["r", "z"]:
        if coord in realistic_wf_map:
            realistic_wf_map[coord] = (
                np.asarray(realistic_wf_map[coord]) * 1000
            ).astype(np.float64)

    # Search for waveform keys (general - I can have any angle in the ideal map)
    keys_to_convolve = [k for k in realistic_wf_map if "waveform" in k]

    for key in keys_to_convolve:
        logger.info("Processing %s...", key)

        # Extract and prepare data
        ideal_wfs = realistic_wf_map[key]
        ideal_wfs_arr = np.asarray(
            ideal_wfs.nda if hasattr(ideal_wfs, "nda") else ideal_wfs
        )
        ideal_wfs_arr = np.nan_to_num(ideal_wfs_arr, nan=0.0)

        original_shape = ideal_wfs_arr.shape
        wfs_flat = ideal_wfs_arr.reshape(-1, original_shape[-1])

        # POST-PROCESSING CHAIN

        # Convolution electronics
        convolved_ak = convolve_waveforms(wfs_flat, rf_kernel)
        wfs_convolved = ak.to_numpy(convolved_ak)

        # Derivative (Charge -> Current)
        wfs_current = to_current(wfs_convolved, dt)

        # Calculate drift time from Amax
        ######################### Maybe it's better to calculate this outside with dedicated function??
        ######################### Should this be calculated after MWA? No noise in the simulations, so this should be accurate?
        current_peak_indices = np.argmax(wfs_current, axis=1)
        drift_indices = current_peak_indices - kernel_delay_idx
        drift_times_flat = drift_indices * dt

        # Add drift time to the realistic map
        dt_key = key.replace("waveform", "drift_time")
        realistic_wf_map[dt_key] = drift_times_flat.reshape(original_shape[:-1]).astype(
            np.float64
        )

        # Moving Window Average
        wfs_mwa = apply_mwa(wfs_current)

        # Align current waveforms to Amax
        wfs_aligned, _ = shift_array(
            wfs_mwa, alignment_idx, nsamples_output_current_wfs
        )

        # Reshape and save
        new_length = wfs_aligned.shape[-1]
        new_shape = (*original_shape[:-1], new_length)
        realistic_wf_map[key] = np.ascontiguousarray(
            wfs_aligned.reshape(new_shape), dtype=np.float64
        )

    return realistic_wf_map


def write_realistic_struct_to_lh5(
    realistic_dict: dict, original_struct: Struct, output_filename: str, group_name: str
) -> Struct:
    """
    Writes realistic current map back to LH5 file as a Struct.
    Output should be compatible with reboost.

    Parameters
    ----------
    realistic_dict : dict
        Dictionary containing processed waveforms, drift times, and
        updated metadata (r, z, t0).
    original_struct : lgdo.Struct
        Original LGDO Struct from the ideal waveform map, used to
        inherit base attributes.
    output_filename : str
        The path to the .lh5 file to be created.
    group_name : str
        The name of the group/Struct inside the LH5 file.

    Returns
    -------
    out_struct : lgdo.Struct
        The final LGDO Struct object that was written to disk.
    """

    struct_data = {}

    for key, val in realistic_dict.items():
        if isinstance(val, np.ndarray):
            # Copy original attributes if they exist
            attrs = (
                copy.deepcopy(original_struct[key].attrs)
                if key in original_struct
                else {}
            )

            # Unit overrides - needed for reboost compatibility
            if "waveform" in key:
                attrs["units"] = "counts"

            elif "drift_time" in key:
                attrs["units"] = "ns"

            elif key in ["r", "z"]:
                attrs["units"] = (
                    "mm"  # m->mm conversion done in realistic map construction
                )

            struct_data[key] = Array(nda=val, attrs=attrs)

        elif hasattr(val, "attrs"):
            struct_data[key] = val

        # New Scalar
        else:
            attrs = {}
            if key == "t0":
                attrs = {"units": "ns", "datatype": "real"}
            struct_data[key] = Scalar(value=val, attrs=attrs)

    # Copy original attributes to preserve simulation metadata
    new_attrs = copy.deepcopy(original_struct.attrs)

    if "datatype" in new_attrs:
        del new_attrs["datatype"]

    # Create LGDO Struct
    out_struct = Struct(obj_dict=struct_data, attrs=new_attrs)

    # Write to disk
    output_path = Path(output_filename)  # Convert string to a Path object

    if output_path.exists():
        logger.info("Overwriting existing file: %s", output_filename)
        output_path.unlink()  # This replaces os.remove()

    lh5.write(obj=out_struct, name=group_name, lh5_file=output_filename, wo_mode="w")

    logger.info("Successfully wrote group '%s' to %s", group_name, output_filename)
    return out_struct
