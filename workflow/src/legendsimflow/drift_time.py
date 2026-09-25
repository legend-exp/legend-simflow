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
"""Compare the simulated and data HPGe drift-time distributions.

The comparison uses the position of the first peak of the distribution and one
of its high percentiles.
"""

from __future__ import annotations

from collections.abc import Mapping

import awkward as ak
import hist
import lh5
import numpy as np
from lgdo import Struct
from matplotlib import pyplot as plt
from numpy.typing import ArrayLike, NDArray
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

from legendsimflow.impurity_tuning import get_grid_value


def get_simulated_drift_times(
    dt_files: list,
    detector: str,
    simid_mapping: Mapping,
    data_stats: Mapping,
    ranges: tuple = (1500, 2500),
) -> tuple[dict[str, ak.Array], Struct]:
    """Get the drift times from the MC for the specified detector and simid mapping.

    This returns the drift times for each run as a dictionary of ak.Arrays and the
    grid info. Only hits with energy in `ranges` are selected. A weight is stored
    based on the number of events in data and simulations.

    Parameters
    ----------
    dt_files
        List of simulation files containing drift times.
    detector
        The detector to read data for.
    simid_mapping
        The mapping from runs to simids (see {func}`get_simid_mapping`).
    data_stats
        Number of events per run in data spectrum, for normalisation.
    ranges
        Range to select drift times.
    """
    drift_time_mc = {}
    for run in simid_mapping:
        files = [file for file in dt_files if simid_mapping[run] in file]

        if len(files) != 1:
            msg = (
                f"Only one drift time file should be present per simid not {len(files)}"
            )
            raise RuntimeError(msg)

        dt_struct = lh5.read(detector, files)

        out = {}
        out["energy"] = dt_struct.energy.view_as("ak")
        psl_scan = {}

        for slope, depv_dict in dt_struct.psl_scan.items():
            psl_scan[slope] = {}
            for dep, dts in depv_dict.items():
                drift_time = dts.drift_time.view_as("ak")
                drift_time = drift_time[
                    (out["energy"] < ranges[1]) & (out["energy"] > ranges[0])
                ]
                psl_scan[slope][dep] = drift_time

        out["energy"] = out["energy"][
            (out["energy"] < ranges[1]) & (out["energy"] > ranges[0])
        ]

        # store a weight for each event based on the number of events in data and simulation
        n_data = data_stats[run]
        n_mc = len(out["energy"])

        if n_mc == 0:
            n_mc = 1

        weights = np.full(n_mc, n_data / n_mc)
        out["weights"] = weights
        out["psl_scan"] = psl_scan

        grid_info = dt_struct.grid_info
        grid_out = {}

        # convert to plain floats
        for a, f in grid_info.items():
            grid_out[a] = float(f.view_as())

        drift_time_mc[run] = ak.Array(out)

    return drift_time_mc, grid_out


def get_simulated_drift_time_obs(
    mc: dict[str, ak.Array],
    grid_info: dict,
    **dt_kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Get the drift time observables from the MC, looping over the grid of parameters.

    Return arrays of:
    - the depletion voltage,
    - the slope,
    - the first drift time observable,
    - the second drift time observable.
    for each combination of slope and depletion voltage in the grid.

    Parameters
    ----------
    mc
        Dictionary of ak.Arrays containing the drift times for each run.
    grid_info
        Information on the parameter grid.
    dt_kwargs
        Keyword arguments to pass to {func}`drift_time_observables`.
    """
    depv = []
    obs1 = []
    obs2 = []
    slopes = []

    # get the grid
    psl_scan = next(iter(mc.values())).psl_scan

    for slope_str in psl_scan.fields:
        slope_idx = int(slope_str.split("_")[-1])
        slope = get_grid_value(slope_idx, grid_info)

        for dep_str in psl_scan[slope_str].fields:
            dep_idx = int(dep_str.split("_")[-1])
            dep = get_grid_value(dep_idx, grid_info, name="dep")

            # get drift times and weights for all runs and concatenate them
            dt = ak.concatenate([mc[run].psl_scan[slope_str][dep_str] for run in mc])
            weights = ak.concatenate([mc[run].weights for run in mc])

            obs = drift_time_observables(dt, weights=weights, **dt_kwargs)[0]

            depv.append(dep)
            slopes.append(slope)
            obs1.append(obs[0])
            obs2.append(obs[1] - obs[0])

    depv = np.array(depv)
    slopes = np.array(slopes)

    obs1 = np.array(obs1)
    obs2 = np.array(obs2)

    return depv, slopes, obs1, obs2


def get_drift_time_chi2(
    data_obs: tuple[np.ndarray, np.ndarray],
    mc_obs: tuple[np.ndarray, np.ndarray],
    weight: float,
) -> np.ndarray:
    """Calculate the chi2 between the data and MC drift time observables."""
    data_obs1, data_obs2 = data_obs
    mc_obs1, mc_obs2 = mc_obs

    return ((data_obs1 - mc_obs1) ** 2 + (data_obs2 - mc_obs2) ** 2) / weight**2


def get_data_drift_times(
    data: dict[str, ak.Array], det: str, ranges=(1500, 2500)
) -> dict[str, NDArray]:
    """Extract the drift time from the data for the specified detector and energy range.

    Returns the drift times as a dictionary keyed by the runid with an array of
    drift times with energy inside `ranges`.

    Parameters
    ----------
    data
        evt tier data per run.
    det
        detector to extract drift time for.
    ranges
        energy range to select drift times.
    """
    energy = {
        run: ak.flatten(d.geds.energy[d.geds.detector_name == det])
        for run, d in data.items()
    }
    drift_time = {
        run: ak.flatten(
            d.geds.psd.low_aoe.time[d.geds.detector_name == det] - d.spms.event_t0
        )
        for run, d in data.items()
    }
    dts = {
        run: drift_time[run][(energy[run] < ranges[1]) & (energy[run] > ranges[0])]
        for run in energy
    }

    return {run: dt[~np.isnan(dt)] for run, dt in dts.items()}


def get_data_drift_time_obs(dts, **kwargs):
    """Get the drift time observables from the data for the specified detector."""
    obs, hist, edges = drift_time_observables(remove_outliers(dts), **kwargs)

    return (obs[0], obs[1] - obs[0]), hist, edges


def plot_drift_time_obs(dts, obs, weights, edges):
    """Plot the drift time observables from the data for the specified detector."""
    h = hist.new.Reg(200, 0, 3200).Double().fill(dts)
    fig, ax = plt.subplots(figsize=(6, 4))
    h.plot(yerr=False)

    h2 = hist.Hist(hist.axis.Variable(edges))
    h2[...] = 16 * weights
    h2.plot(yerr=False)

    ax.set_xlabel("Drift time [ns]")
    ax.set_ylabel("Counts")

    ax.axvline(obs[0], label="Mode", linestyle="--", color="black")
    ax.axvline(obs[1] + obs[0], label="Q-90", linestyle="--", color="red")
    ax.legend()
    return fig


def drift_time_observables(
    drift_times: ArrayLike,
    *,
    weights: ArrayLike | None = None,
    percentile: float = 90,
    smoothing: float = 50,
    peak_threshold: float = 0.5,
) -> tuple[NDArray, NDArray, NDArray]:
    """First peak and a high percentile of a drift-time distribution.

    Both are read from the histogram of the drift times in 1 ns bins. The first
    peak is the lowest drift-time peak higher than `peak_threshold` times the
    global maximum of the histogram smoothed with a Gaussian of standard
    deviation `smoothing`. The percentile is the lower edge of the bin where the
    cumulative histogram reaches it.



    Parameters
    ----------
    drift_times
        drift time of each event in ns. NaNs are ignored.
    weights
        weight of each event, all equal by default.
    percentile
        percentile of the distribution, in percent.
    smoothing
        standard deviation of the Gaussian in ns.
    peak_threshold
        minimum peak height, as a fraction of the global maximum. ``1`` selects
        the global maximum.

    Returns
    -------
    obs
        ``[peak, percentile]`` in ns.
    density
        smoothed histogram of the drift times.
    edges
        bin edges of the histogram.
    """
    x = np.asarray(drift_times, dtype=float)

    finite = np.isfinite(x)
    x = x[finite]
    w = None if weights is None else np.asarray(weights, dtype=float)[finite]

    # padded so that a peak at the edge of the data is found too
    edges = np.arange(x.min() - smoothing, x.max() + smoothing + 1)
    hist = np.histogram(x, edges, weights=w)[0].astype(float)
    density = gaussian_filter1d(hist, smoothing, mode="constant")
    peak = edges[find_peaks(density, height=peak_threshold * density.max())[0][0]]

    cdf = np.cumsum(hist)
    q = edges[np.searchsorted(cdf, percentile / 100 * cdf[-1])]
    return np.array([peak + 0.5, q]), density, edges


def drift_time_cost(
    data_drift_times: ArrayLike,
    sim_drift_times: ArrayLike,
    sim_weights: ArrayLike | None = None,
    **obs_kwargs,
) -> float:
    """Distance between the data and simulated drift-time distributions.

    Sum of the squared differences (ns²) between the
    :func:`drift_time_observables` of the simulation and of the data.

    Parameters
    ----------
    data_drift_times
        drift time of each data event in ns.
    sim_drift_times
        drift time of each simulated event in ns.
    sim_weights
        weight of each simulated event, all equal by default.
    **obs_kwargs
        passed to :func:`drift_time_observables`.
    """
    diff = drift_time_observables(sim_drift_times, weights=sim_weights, **obs_kwargs)[0]
    diff -= drift_time_observables(data_drift_times, **obs_kwargs)[0]
    return float(np.sum(diff**2))


def remove_outliers(values: ArrayLike, percentile: float = 99) -> NDArray:
    """Values up to the `percentile` (in percent) of the sample; NaNs are dropped."""
    x = np.asarray(values)
    return x[x <= np.nanpercentile(x, percentile)]
