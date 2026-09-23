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


import argparse
import logging
from collections.abc import Mapping
from pathlib import Path

import awkward as ak
import dbetto
import hist
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
import lh5
import matplotlib.pyplot as plt
import numpy as np
from lgdo import Struct
from matplotlib.backends.backend_pdf import PdfPages
from numpy.typing import ArrayLike
from scipy.interpolate import griddata
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import drift_time, utils
from legendsimflow.metadata import get_simconfig
from legendsimflow.plot import decorate
from legendsimflow.scripts import log_script_invocation

DEFAULT_SETTINGS = {
    "drift_time_weight": 50,  # ns
    "wf_weight": 0.5,  # arb
    "dt_kwargs": {"percentile": 90, "smoothing": 50, "peak_threshold": 0.25},
}


def get_run_mapping(simconfig: Mapping, runs) -> dict[str, str]:
    """Get a mapping from run to simid from the simconfig."""
    out = {}
    for simid, info in simconfig.items():
        if len(info["runlist"]) != 1:
            msg = f"Only one run per simid is supported. Found {len(info.runlist)} runs for simid {simid}."
            raise ValueError(msg)

        if any(run in info.runlist[0] for run in runs):
            out[info.runlist[0]] = simid
    return out


def _get_grid_value(idx: int, grid_info: dict, name="slope") -> float:
    """Get the slope value from the grid info given the name."""
    return float(grid_info[f"{name}_min"] + grid_info[f"{name}_step"] * idx)


def get_drift_times_mc(dt_files, det, simid_mapping, run_norms):
    drift_time_mc = []
    weights = []

    for run in simid_mapping:
        files = [file for file in dt_files if simid_mapping[run] in file]
        if len(files) != 1:
            msg = (
                f"Only one drift time file should be present per simid not {len(files)}"
            )
            raise RuntimeError(msg)

        drift_times = lh5.read(det, files)
        weight = ak.full_like(drift_times.energy.view_as("ak"), run_norms[run])

        weights.append(weight)
        drift_time_mc.append(drift_times)
    return weights, drift_time_mc


def get_drift_time_obs_mc(
    mc: Struct, grid_info: dict, weights: ArrayLike, ranges=(1500, 2500), **dt_kwargs
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Get the drift time observables from the MC."""
    energy = [m.energy.view_as("ak") for m in mc]
    weights = ak.concatenate(
        [
            w[(e > ranges[0]) & (e < ranges[1])]
            for w, e in zip(weights, energy, strict=True)
        ]
    )
    depv = []
    obs1 = []
    obs2 = []
    slopes = []

    for slope_str in mc[0].psl_scan:
        slope_idx = int(slope_str.split("_")[-1])
        slope = _get_grid_value(slope_idx, grid_info)

        for dep_str in mc[0].psl_scan[slope_str]:
            dep_idx = int(dep_str.split("_")[-1])
            dep = _get_grid_value(dep_idx, grid_info, name="dep")

            dt = ak.concatenate(
                [
                    out.psl_scan[slope_str][dep_str].drift_time.view_as("ak")[
                        (e > ranges[0]) & (e < ranges[1])
                    ]
                    for e, out in zip(energy, mc, strict=True)
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


def get_dt_chi2(
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
        files = list(Path(path_data).glob(f"*-{run}-*.lh5"))

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


def get_drift_time(data: dict[str, ak.Array], det: str, ranges=(1500, 2500)):
    """Extract the drift time from the data for the specified detector and energy range."""
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
    dts = np.concatenate(
        [
            drift_time[run][(energy[run] < ranges[1]) & (energy[run] > ranges[0])]
            for run in energy
        ]
    )

    return dts[~np.isnan(dts)]


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


def plot_surface(x, y, z, name, det, vrange, levels, method="nearest"):
    """Plot the cost function as a function of depletion voltage and slope."""
    xi = np.linspace(x.min(), x.max(), 500)
    yi = np.linspace(y.min(), y.max(), 500)
    X, Y = np.meshgrid(xi, yi)

    Z = griddata((x, y), z, (X, Y), method=method)

    fig, ax = plt.subplots()

    cmap = plt.colormaps["RdYlBu_r"].copy()
    cmap.set_over("grey")

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
    plt.show()

    return fig, ax, xm, ym, zm


@snakemake_compatible(
    mapping={
        "elecmod": "input.elecmod",
        "drift_time": "input.drift_time",
        "data_files": "input.data_files",
        "pars_file": "output.pars_file",
        "plot_file": "output.plot_file",
        "settings": "input.settings",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract the HPGe electronics model for a LEGEND run."
    )

    parser.add_argument(
        "--elecmod",
        help="input YAML file for the electronics model",
    )
    parser.add_argument(
        "--drift-time-files",
        required=True,
        nargs="+",
        help="input LH5 file for the drift times.",
    )
    parser.add_argument(
        "--pars-file",
        required=True,
        help="output YAML file for the electronics model parameters",
    )
    parser.add_argument(
        "--data-path",
        required=True,
        help="input path for evt tier data.",
    )
    parser.add_argument(
        "--runs",
        required=True,
        nargs="+",
        help="list of runs to process.",
    )
    parser.add_argument(
        "--run-norms",
        required=True,
        help="input YAML file for the run normalisation (e.g. livetime).",
    )

    parser.add_argument("--log-file", default=None, help="log file")
    parser.add_argument(
        "--simflow-config",
        "--config",
        dest="simflow_config",
        required=True,
        help="simflow config YAML path",
    )
    parser.add_argument(
        "--settings",
        type=str,
        required=False,
        default=None,
        help="Path to YAML file with settings for the fit (e.g. initial values, limits, comparison window); ",
    )
    parser.add_argument(
        "--plot-file",
        type=str,
        required=False,
        default=None,
        help="File name for diagnostic plots.",
    )

    args = parser.parse_args()
    config = utils.init_simflow_context(args.simflow_config, workflow=None).config
    metadata = config.metadata

    if args.log_file is not None:
        log = ldfs.utils.build_log(metadata.simprod.config.logging, args.log_file)
    else:
        logging.basicConfig(
            level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
        )
        log = logging.getLogger(__name__)

    settings = (
        dbetto.AttrsDict(dbetto.utils.load_dict(args.settings))
        if args.settings is not None
        else dbetto.AttrsDict(DEFAULT_SETTINGS)
    )
    log_script_invocation(log, "extract-hpge-impurity-model", parser, args)

    # 1. load data
    msg = f"... loading data from runs {args.runs} and {args.data_path}"
    log.info(msg)
    data = read_data(args.data_path, args.runs)

    # 2. get run norms (e.g. from livetime)
    log.info("... loading run norms")
    run_norms = dbetto.utils.load_dict(args.run_norms)

    simid_mapping = get_run_mapping(get_simconfig(config, "hit", simid=None), args.runs)

    out = {}
    dets = lh5.ls(args.drift_time_files[0])

    with PdfPages(args.plot_file) as pdf:
        for det in dets:
            msg = f"... processing {det}"
            log.info(msg)
            grid_info = {
                a: f.view_as()
                for a, f in lh5.read(
                    f"{det}/grid_info", args.drift_time_files[0]
                ).items()
            }
            # 3. get data observables
            dts = get_drift_time(data, det)
            data_dt_obs, weights, edges = get_drift_time_obs(dts, **settings.dt_kwargs)

            fig = plot_drift_time_obs(dts, data_dt_obs, weights, edges)
            decorate(fig)
            pdf.savefig()
            plt.close(fig)

            log.info("... found data observables (%f, %f)", *data_dt_obs)

            # 4. get mc observables
            weights, dt_mc = get_drift_times_mc(
                args.drift_time_files, det, simid_mapping, run_norms
            )

            depv, slope, dt_obs1, dt_obs2 = get_drift_time_obs_mc(
                dt_mc, grid_info, weights, **settings.dt_kwargs
            )
            dt_chi2 = get_dt_chi2(
                data_dt_obs, (dt_obs1, dt_obs2), settings.drift_time_weight
            )

            fig, _, best_dep, best_slope, best_cost = plot_surface(
                depv,
                slope,
                dt_chi2,
                r"$\chi^2$",
                det,
                vrange=(0, 10),
                levels=[2, 5, 10],
                method="nearest",
            )

            decorate(fig)
            pdf.savefig()
            plt.close(fig)

            msg = f"For {det} found minimum Vdep = {best_dep:1f}, slope {best_slope:.1f} with chi2 {best_cost:.2f}"
            log.info(msg)

            # 5. extract electronics model parameters

            # elecmod = dbetto.utils.load_dict(args.elecmod)

            # wf_chi2 = get_wf_chi2(elecmod,settings.wf_weight)

            # find the best fit
            # best_slope, best_dep  = fit_impurities(det,wf_chi2,dt_chi2,pdf)

            out[det] = {"slope": best_slope, "depletion_voltage": best_dep}

    dbetto.utils.write_dict(out, args.pars_file)


if __name__ == "__main__":
    main()
