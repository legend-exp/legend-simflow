"""
Suggested changes:
- Read convolution parameters from some sort of metadata file instead of command line
"""


import os
import argparse
import numpy as np
import awkward as ak
import copy
from scipy.signal import convolve, fftconvolve
import pygama
from lgdo import lh5, Struct, Array, Scalar
from dspeed.processors import moving_window_multi


ALIGNMENT_IDX = 2000  # Index to align current waveforms to Amax
######################### Want to put 1000 as default???

NSAMPLES_OUTPUT_CURRENT_WFS = 4001  # Final length of the realistic current waveforms in the map


def create_response_kernel(dt, mu_bandwidth, sigma_bandwidth, tau_rc, gaussian_only=False):
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
    # Convert to samples
    mu_samples = mu_bandwidth / dt
    sigma_samples = sigma_bandwidth / dt
    tau_samples = tau_rc / dt
    
    # Define the sample range
    x = np.arange(-100, 500)
    
    # Compute Gaussian response (digitizer bandwidth)
    rf_digi = np.exp(-0.5 * ((x - mu_samples) / sigma_samples)**2)
    rf_digi /= np.sum(rf_digi)
    
    if gaussian_only:
        return rf_digi
    
    # Compute causal exponential decay (preamplifier response) - response is zero for x < 0
    rf_preamp = np.zeros_like(x, dtype=float)
    rf_preamp[x >= 0] = np.exp(-x[x >= 0] / tau_samples)
    rf_preamp /= np.sum(rf_preamp)

    # Convolve preamp and digitizer responses - 'full' mode results in a length of (len(rf_preamp) + len(rf_digi) - 1)
    rf = convolve(rf_preamp, rf_digi, mode='full')
    
    return rf


def convolve_waveforms(wf_array, rf_kernel):
    """
    Convolve an array of waveforms with the response kernel
    
    Parameters
    ----------
    wf_array : ak.Array
        Array of waveforms (each waveform is a 1D array)
    rf_kernel : np.ndarray
        The response kernel from create_response_kernel()
    
    Returns
    -------
    convolved_wfs : ak.Array
        Array of convolved waveforms with same shape as input
    """
    convolved_list = []
    
    for wf in wf_array:
        # Convert to numpy array
        padded_wf = np.asarray(wf, dtype=float).flatten()
        
        # Convolve with kernel using 'full' mode
        wf_conv = convolve(padded_wf, rf_kernel, mode='full')

        # Slicing - correct???
        # Slice to original length to maintain timing relative to start.
        # This preserves the 'warm-up' phase to emulate baseline
        wf_conv = wf_conv[:len(padded_wf)]  
        
        convolved_list.append(wf_conv)
    
    return ak.Array(convolved_list)


def convolve_waveforms_fast(wf_array, rf_kernel, batch_size=50000):
    """
    Vectorized convolution using FFT with batching to save memory.
    
    Parameters
    ----------
    wf_array : ak.Array
        Array of waveforms (must be same length!!!)
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
        # Convert batch to 2D numpy matrix
        batch = wf_array[i : i + batch_size]
        wf_matrix = ak.to_numpy(batch)
        
        # FFT convolution 
        wf_conv = fftconvolve(wf_matrix, kernel_2d, mode='full')
        
        # Slicing - correct???
        # Slice to original length to maintain timing relative to start.
        # This preserves the initial points to emulate baseline
        wf_conv_sliced = wf_conv[:, :wf_matrix.shape[1]]
        
        convolved_results.append(wf_conv_sliced)

    # Concatenate batches and return as Awkward Array
    return ak.concatenate(convolved_results)



def to_current(wf_array, dt):
    """
    Computes the derivative of each charge waveform 

    Parameters:
    -----------
    wf_simulated : numpy.ndarray
        2D array of waveforms
    dt : float
        Time step between samples (in ns!!)
    
    Returns:
    -----------
    derivatives : numpy.ndarray
        2D array containing the derivatives - each has same shape as input!!
    """
    ########################## Normalization to sampling rate ?!? 
    derivatives = np.diff(wf_array, axis=1, prepend=0) / dt
    return derivatives


def apply_mwa(wf_array):
    """
    Apply moving window average to each current waveform in the array.

        Parameters
    ----------
    wf_array : np.ndarray
        2D Array of current waveforms [n_events, n_samples]
    
    Returns
    -------
    mwa_currents : np.ndarray
        Array of filtered waveforms - same shape as input!
    """

    mwa_currents = np.zeros_like(wf_array, dtype=float) # Risk to saturate memory?? Should check this??

    ########## Should vectorize?? ########## 
    # If run on map only, this is fast enough
    for i in range(len(wf_array)):
        wf_in = np.ascontiguousarray(wf_array[i])
        wf_out = mwa_currents[i]
        ########## Make sure these are the best parameters!!
        moving_window_multi(wf_in, 48, 3, 0, wf_out)
        
    return mwa_currents



def shift_array(wf_input, center_idx=ALIGNMENT_IDX, output_len=NSAMPLES_OUTPUT_CURRENT_WFS):
    """
    Aligns an array of waveforms by shifting their maximum to a fixed index.
    No normalization is performed; raw amplitudes are preserved.

    Parameters:
    -----------
    wf_input : ak.Array or np.ndarray
        Input array of current waveforms 
    center_idx : int
        The index in the output array where the peak will be placed 
    output_len : int
        The total length of the resulting aligned current waveforms 
        
    Returns:
    --------
    shifted_wfs : np.ndarray
        2D array of shifted waveforms.
    peak_indices : np.ndarray
        1D array containing the original peak index for each current waveform
        ######################### It is NOT the original drift time as current waveform inherits baseline from convolution!!!!!!
    """
    
    # Convert to Numpy for matrix operations
    if isinstance(wf_input, ak.Array):
        wfs = ak.to_numpy(wf_input)
    else:
        wfs = wf_input


    # Output Array - Indeces not covered by the shifted signal are automatically padded with 0.0
    n_wfs, n_samples = wfs.shape
    shifted_wfs = np.zeros((n_wfs, output_len), dtype=float)
    
    # Find Amax
    peak_indices = np.argmax(wfs, axis=1)
    
    # Shift each waveform
    for i in range(n_wfs):
        peak = peak_indices[i]
        shift = center_idx - peak
        
        # Bounds to prevent IndexErrors
        start_src = max(0, -shift)
        end_src   = min(n_samples, output_len - shift)
        
        if start_src < end_src:
            start_dst = start_src + shift
            end_dst   = end_src   + shift
        
            shifted_wfs[i, start_dst:end_dst] = wfs[i, start_src:end_src]

    return shifted_wfs, peak_indices



def generate_realistic_map(ideal_wf_map_obj, rf_kernel, Amax_alignment_idx=ALIGNMENT_IDX):
    """
    Applies the waveform post-processing chain to generate a realistic waveform map starting from an ideal one:
    1. Converts coordinates (m to mm) ######################### Correct? Ideal wf map is in meters while reboost wants mm????
    2. Convolves with system response
    3. Aligns by Peak Time 
    4. Calculates compensated Drift Time
    """
    
  
    realistic_wf_map = dict(ideal_wf_map_obj)
    dt = ideal_wf_map_obj['dt'].value if hasattr(ideal_wf_map_obj['dt'], 'value') else ideal_wf_map_obj['dt']
    realistic_wf_map['t0'] = -1.0 * Amax_alignment_idx * dt # Set the global t0 relative to alignment index
    
    # Delay of the response kernel
    kernel_delay_idx = np.argmax(rf_kernel)

    # COORDINATE FIX: Convert m to mm ######################### Correct? 
    for coord in ['r', 'z']:
        if coord in realistic_wf_map:
            realistic_wf_map[coord] = (np.asarray(realistic_wf_map[coord]) * 1000).astype(np.float64)

    # Search for waveform keys (general - I can have any angle in the ideal map)
    keys_to_convolve = [k for k in realistic_wf_map.keys() if "waveform" in k]

    for key in keys_to_convolve:
        print(f"Processing {key}...")
        
        # Extract and prepare data
        ideal_wfs = realistic_wf_map[key]
        ideal_wfs_arr = np.asarray(ideal_wfs.nda if hasattr(ideal_wfs, 'nda') else ideal_wfs)
        ideal_wfs_arr = np.nan_to_num(ideal_wfs_arr, nan=0.0)

        original_shape = ideal_wfs_arr.shape
        wfs_flat = ideal_wfs_arr.reshape(-1, original_shape[-1])
        
        # -----------  POST-PROCESSING CHAIN
        
        # Convolution electronics
        convolved_ak = convolve_waveforms(wfs_flat, rf_kernel)
        wfs_convolved = ak.to_numpy(convolved_ak)
        
        # Derivative (Charge -> Current)
        wfs_current = to_current(wfs_convolved, dt)
        
        # Calculate drift time from Amax (Compensated for Convolution Delay)
        ######################### Maybe it's better to calculate this outside with dedicated function??
        ######################### Should this be calculated after MWA? No noise in the simulations, so this should be accurate?
        current_peak_indices = np.argmax(wfs_current, axis=1)
        drift_indices = current_peak_indices - kernel_delay_idx
        drift_times_flat = drift_indices * dt

        # Add drfift time to the realistic map
        dt_key = key.replace("waveform", "drift_time")
        realistic_wf_map[dt_key] = drift_times_flat.reshape(original_shape[:-1]).astype(np.float64)

        # Moving Window Average 
        wfs_mwa_ak = apply_mwa(wfs_current)
        wfs_mwa = ak.to_numpy(wfs_mwa_ak) 
        
        # Align current waveforms to Amax
        wfs_aligned, _ = shift_array(wfs_mwa, center_idx=Amax_alignment_idx)

        # Reshape and save
        new_length = wfs_aligned.shape[-1]
        new_shape = original_shape[:-1] + (new_length,)
        realistic_wf_map[key] = np.ascontiguousarray(wfs_aligned.reshape(new_shape), dtype=np.float64)
        
    return realistic_wf_map




def write_realistic_struct_to_lh5(realistic_dict, original_struct, output_filename, group_name):
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
            attrs = copy.deepcopy(original_struct[key].attrs) if key in original_struct else {}
            
            # Unit overrides ####### Needed for reboost compatibility
            if "waveform" in key:
                attrs['units'] = 'counts'
            
            elif "drift_time" in key:
                attrs['units'] = 'ns'
            
            elif key in ['r', 'z']:
                attrs['units'] = 'mm' # m->mm conversion done in realistic map construction
            
            struct_data[key] = Array(nda=val, attrs=attrs)
            
        elif hasattr(val, 'attrs'): 
            struct_data[key] = val
            
        # New Scalar
        else:
            attrs = {}
            if key == 't0':
                attrs = {'units': 'ns', 'datatype': 'real'}
            struct_data[key] = Scalar(value=val, attrs=attrs)


    # Copy original attributes to preserve simulation metadata
    new_attrs = copy.deepcopy(original_struct.attrs)
    
    if 'datatype' in new_attrs:
        del new_attrs['datatype']

    # Create LGDO Struct
    out_struct = Struct(obj_dict=struct_data, attrs=new_attrs)

    # Write to disk
    if os.path.exists(output_filename):
        print(f"Overwriting existing file: {output_filename}")
        os.remove(output_filename)
        
    lh5.write(obj=out_struct, name=group_name, lh5_file=output_filename, wo_mode='w')
    
    print(f"Successfully wrote group '{group_name}' to {output_filename}")
    return out_struct





def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--detector', required=True, help="Detector name (LH5 group)")
    parser.add_argument('--sigma-conv', type=float, required=True, help="Sigma of the gaussian component of the convolution kernel in ns")
    parser.add_argument('--tau-conv', type=float, required=True, help="Tau of the exponential component of the convolution kernel in ns")
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--output-file', required=True)
    args = parser.parse_args()

    # 1. Load data
    ideal_map_obj = lh5.read(f"{args.detector}/", args.input_file)
    dt = ideal_map_obj['dt'].value
    
    # 2. Setup Physics Kernel (mu=0, sigma=15ns, tau=70ns)
    rf_kernel = create_response_kernel(dt, mu_bandwidth=0, sigma_bandwidth=args.sigma_conv, tau_rc=args.tau_conv)

    # 3. Process
    realistic_dict = generate_realistic_map(ideal_map_obj, rf_kernel)

    # 4. Write output with units
    write_realistic_struct_to_lh5(
        realistic_dict=realistic_dict,
        original_struct=ideal_map_obj,
        output_filename=args.output_file,
        group_name=args.detector
    )
    print(f"Realistic library created successfully: {args.output_file}")

if __name__ == "__main__":
    main()