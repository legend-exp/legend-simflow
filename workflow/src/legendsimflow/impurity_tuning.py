# ruff: noqa: I002

# Copyright (C) 2026 Toby Dixon <toby.dixon.23@ucl.ac.uk>
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

from collections.abc import Mapping
from pathlib import Path

import awkward as ak
import hist
import lh5
import matplotlib.pyplot as plt
import numpy as np
from lgdo import Struct
from numpy.typing import NDArray
from scipy.interpolate import griddata

from legendsimflow import drift_time

DEFAULT_SETTINGS = {
    "drift_time_weight": 50,  # ns
    "wf_weight": 0.5,  # arb
    "dt_kwargs": {"percentile": 90, "smoothing": 50, "peak_threshold": 0.25},
}


def get_simid_mapping(simconfig: Mapping, runs: list[str]) -> dict[str, str]:
    """Get a mapping from run to simid from the simconfig.

    Only simids with a single run per simid are supported, and only
    runs present in the `runs` list are used.

    Returns a dictionary mapping runs to simids.

    Parameters
    ----------
    simconfig
        metadata on simflow configuration including runlists for each simid.
    runs
        lost of runids to use.
    """
    out = {}
    for simid, info in simconfig.items():
        if not any(run in info.runlist[0] for run in runs):
            continue

        if len(info["runlist"]) != 1:
            msg = f"Only one run per simid is supported. Found {len(info.runlist)} runs for simid {simid}."
            raise ValueError(msg)

        out[info.runlist[0]] = simid

    return out


def get_grid_value(idx: int, grid_info: dict, name="slope") -> float:
    """Get the slope value from the grid info given the name.

    Parameters
    ----------
    idx
        The index of the grod to extract
    grid_info
        Information on the grid, must contain `{name}_min` and {name}_step
    name
        The field to extract.
    """
    return float(grid_info[f"{name}_min"] + grid_info[f"{name}_step"] * idx)


def get_simulated_drift_times(
    dt_files: list,
    detector: str,
    simid_mapping: Mapping,
    data_stats: Mapping,
    ranges: tuple = (1500, 2500),
) -> tuple[dict[str, ak.Array],Struct]:
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
        
        for slope, depv_dict in dt_struct.psl_scan.items():
            out[slope ] = {}
            for dep, dts in depv_dict.items():

                drift_time = dts.drift_time.view_as("ak")
                drift_time = drift_time[(out["energy"]<ranges[1]) & (out["energy"]>ranges[0])]
                out[slope][dep] = drift_time

        # store a weight for each event based on the number of events in data and simulation
        n_data = data_stats[run]
        n_mc = len(out["energy"])
        
        if n_mc == 0:
            n_mc = 1

        weights = np.full(n_mc, n_data / n_mc)
        out["weights"] = weights

        grid_info = dt_strict.grid_info
        drift_time_mc[run] = ak.Array(out)

    return drift_time_mc, grid_info

def get_simulated_drift_time_obs(
    mc: Struct,
    grid_info: dict,
    **dt_kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Get the drift time observables from the MC, looping over the grid of parameters."""
    depv = []
    obs1 = []
    obs2 = []
    slopes = []

    energy = [out.energy.view_as("ak") for out in mc.values()]

    # get the correct weights
    weights = [
        w[(e > ranges[0]) & (e < ranges[1])]
        for e, w in zip(energy, weights.values(), strict=True)
    ]
    weights = np.concatenate([w / len(w) for w in weights])
    psl_scan = next(iter(mc.values())).psl_scan

    for slope_str in psl_scan:
        slope_idx = int(slope_str.split("_")[-1])
        slope = get_grid_value(slope_idx, grid_info)

        for dep_str in psl_scan[slope_str]:
            dep_idx = int(dep_str.split("_")[-1])
            dep = get_grid_value(dep_idx, grid_info, name="dep")

            dt = ak.concatenate(
                [
                    out.psl_scan[slope_str][dep_str].drift_time.view_as("ak")[
                        (e > ranges[0]) & (e < ranges[1])
                    ]
                    for e, out in zip(energy, mc.values(), strict=True)
                ]
            )
            obs = drift_time.drift_time_observables(dt, weights=weights, **dt_kwargs)[0]

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


def read_data(path_data: str, runs: list[str]) -> dict[str, ak.Array]:
    """Read the evt tier data from the specified runs and return a dictionary of ak.Arrays."""
    data = {}

    for run in runs:
        files = list(Path(path_data).glob(f"*{run}*.lh5"))

        if len(files) == 0:
            msg = "No data files found!"
            raise RuntimeError(msg)

        data[run] = lh5.read(
            "evt",
            files,
            field_mask=[
                "coincident",
                "geds/energy",
                "geds/detector_name",
                "geds/psd/low_aoe",
                "trigger",
                "spms/event_t0",
            ],
        ).view_as("ak")

        # TODO. some cuts

    return data


def get_drift_time(
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


def get_drift_time_obs(dts, **kwargs):
    """Get the drift time observables from the data for the specified detector."""
    obs, hist, edges = drift_time.drift_time_observables(
        drift_time.remove_outliers(dts), **kwargs
    )

    return (obs[0], obs[1] - obs[0]), hist, edges


def plot_drift_time_obs(dts, obs, weights, edges):
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


def plot_cost_surface(x, y, z, name, det, vrange, levels, method="nearest"):
    """Plot the cost function as a function of depletion voltage and slope."""
    xi = np.linspace(x.min(), x.max(), 500)
    yi = np.linspace(y.min(), y.max(), 500)
    X, Y = np.meshgrid(xi, yi)

    Z = griddata((x, y), z, (X, Y), method=method)

    fig, ax = plt.subplots()

    cmap = plt.colormaps["RdYlBu_r"].copy()
    cmap.with_extremes(over="grey")

    im = ax.pcolormesh(
        X,
        Y,
        Z,
        cmap=cmap,
        vmin=vrange[0],
        vmax=vrange[1],
        shading="auto",
        rasterized=True,
    )

    xm = x[np.argmin(z)]
    ym = y[np.argmin(z)]
    zm = np.min(z)
    ax.scatter([xm], [ym], color="red", s=50)

    fig.colorbar(im, ax=ax, label=name, cmap="cividis")

    # Values at which to draw contours
    cs = ax.contour(
        X,
        Y,
        Z,
        levels=levels,
        colors="black",
    )
    ax.clabel(cs, fmt="%g", fontsize=12)
    ax.set_xlabel("Depletion voltage [V]")
    ax.set_ylabel("Slope [%]")
    ax.set_title(f"{det} cost function (min {np.min(z):.2f})")

    return fig, ax, xm, ym, zm
