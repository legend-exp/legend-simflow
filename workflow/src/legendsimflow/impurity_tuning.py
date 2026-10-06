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

import awkward as ak
import lh5
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from .utils import get_evt_tier_name, lookup_evt_files


def get_grid_value(idx: int, grid_info: dict, name="slope") -> float:
    """Get the slope value from the grid info given the name.

    Parameters
    ----------
    idx
        The index of the grod to extract
    grid_info
        Information on the grid, must contain `{name}_min` and {name}_step.
    name
        The field to extract.
    """
    return float(grid_info[f"{name}_min"] + grid_info[f"{name}_step"] * idx)


def read_evt_data(path_data: str, runs: list[str]) -> dict[str, ak.Array]:
    """Read the evt tier data from the specified runs and return a dictionary of ak.Arrays.

    Parameters
    ----------
    path_data
        Path to the data files.
    runs
        List of runids to read data for.

    Returns
    -------
    data
        Dictionary of ak.Arrays containing the evt data for each run.
    """
    data = {}
    evt_tier_name = get_evt_tier_name(path_data)

    for runid in runs:
        files = lookup_evt_files(path_data, runid, evt_tier_name)

        if len(files) == 0:
            msg = f"no {evt_tier_name} tier files found for {runid} in {path_data}"
            raise RuntimeError(msg)

        evt_data = lh5.read(
            "evt",
            files,
            field_mask=[
                "coincident",
                "geds/energy",
                "geds/detector_name",
                "geds/psd/low_aoe",
                "geds/quality",
                "trigger",
                "spms/event_t0",
                "spms/energy_sum",
            ],
        ).view_as("ak")

        mask = (
            ak.all(evt_data.geds.quality.is_good_channel, axis=-1)
            & (~evt_data.trigger.is_forced)
            & (~evt_data.coincident.puls)
            & (~evt_data.coincident.muon)
            & (~evt_data.coincident.muon_offline)
            & evt_data.geds.quality.is_bb_like
            & (evt_data.spms.energy_sum > 10)
        )
        data[runid] = evt_data[mask]

    return data


def plot_cost_surface(
    x, y, z, name, det, vrange, levels, method="nearest", ftype="drift time"
):
    """Plot the cost function as a function of depletion voltage and slope."""
    xi = np.linspace(x.min(), x.max(), 500)
    yi = np.linspace(y.min(), y.max(), 500)
    X, Y = np.meshgrid(xi, yi)

    Z = griddata((x, y), z, (X, Y), method=method)

    fig, ax = plt.subplots()

    cmap = plt.colormaps["RdYlBu_r"].copy()
    cmap = cmap.with_extremes(over="grey")

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
    ax.set_title(f"{det} {ftype} (min {np.min(z):.2f})")

    return fig, ax, xm, ym, zm


def get_wf_chi2(
    elecmod: Mapping, grid_info: Mapping, wf_scale
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Get the waveform chi2 for the given electronics model.

    Parameters
    ----------
    elecmod
        Electronics model parameters.
    grid_info
        Information on the grid, must contain `dep_min`, `dep_step`, `slope_min`, and `slope_step`.
    wf_scale
        Scale factor for the RMS to chi2.

    Returns
    -------
    dep
        Depletion voltage values.
    slope
        Slope values.
    wf_chi2
        Waveform chi2 values.
    """
    dep = []
    slope = []
    wf_rms = []
    
    for slope_val, slope_dict in elecmod.items():
        for depv_val, depv_dict in slope_dict.items():
            
            dep.append(get_grid_value(float(depv_val.split("_")[-1]), grid_info, name="dep"))
            slope.append(get_grid_value(float(slope_val.split("_")[-1]), grid_info, name="slope"))
            wf_rms.append(depv_dict["rms"] ** 2 / wf_scale**2)

    return dep, slope, wf_rms
