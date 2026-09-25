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
import lh5
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

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

    for run in runs:
        files = list(Path(path_data).glob(f"*{run}*.lh5"))

        if len(files) == 0:
            msg = "No data files found!"
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
        data[run] = evt_data[mask]

    return data


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
