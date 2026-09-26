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

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
import lh5
import matplotlib.pyplot as plt
import numpy as np
import reboost
from matplotlib.backends.backend_pdf import PdfPages
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, utils
from legendsimflow.drift_time import (
    get_data_drift_time_obs,
    get_data_drift_times,
    get_drift_time_chi2,
    get_simulated_drift_time_obs,
    get_simulated_drift_times,
    plot_drift_time_obs,
)
from legendsimflow.impurity_tuning import plot_cost_surface, read_evt_data
from legendsimflow.metadata import get_runlist
from legendsimflow.plot import decorate
from legendsimflow.scripts import log_script_invocation

DEFAULT_SETTINGS = {
    "drift_time_weight": 50,  # ns
    "wf_weight": 0.5,  # arb
    "dt_kwargs": {"percentile": 90, "smoothing": 50, "peak_threshold": 0.25},
    "energy_range": [1500, 2500],
}


@snakemake_compatible(
    mapping={
        "drift_time_files": "input.drift_time",
        "data_path": "params.data_path",
        "pars_file": "output.pars_file",
        "plot_file": "output.plot_file",
        "settings": "input.settings",
        "simids": "params.simids",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    """Extract the HPGe impurity curve model for a LEGEND run.

    - This script reads the data and extracts the drift time observables from the data and MC.
    - It then calculates the chi2 between the data and MC observables and finds the best fit impurity curve scalings.
    - Finally, it writes the best fit parameters to a YAML file and produces diagnostic plots.

    The output is a YAML file with the following format:
    ```yaml
    detector_name:
        slope: <best fit slope>
        depletion_voltage: <best fit depletion voltage>
    ```

    The simulations are required to have a specific format with one simid per run.

    The normalisation is fine determined per run, before combination.

    """
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
        "--simids",
        required=True,
        nargs="+",
        help="simulation IDs of the drift-time files, in the same order.",
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
        required=True,
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

    data_path = nersc.dvs_ro(config, args.data_path)
    drift_time_files = nersc.dvs_ro(config, args.drift_time_files)

    # make a profiler to track performance of the script
    perf_block, print_perf, _ = reboost.make_profiler()

    # the settings file is shared with the drift-time scan, so it may hold only
    # part of the fit settings
    settings = dbetto.AttrsDict(DEFAULT_SETTINGS)
    if args.settings is not None:
        settings = dbetto.AttrsDict(
            DEFAULT_SETTINGS | dbetto.utils.load_dict(args.settings)
        )
    log_script_invocation(log, "extract-hpge-impurity-model", parser, args)

    # each simulation holds a single run, which its drift times are compared to
    run_files = {}
    for simid, file in zip(args.simids, drift_time_files, strict=True):
        runs = get_runlist(config, simid)
        if len(runs) != 1:
            msg = f"simid {simid} must have a single run in its runlist, found {runs}"
            raise ValueError(msg)
        run_files[runs[0]] = file

    # 1. load data
    msg = f"... loading data from runs {list(run_files)} and {args.data_path}"
    log.info(msg)

    with perf_block("read_evt_data()"):
        data = read_evt_data(data_path, list(run_files))

    out = {}
    dets = lh5.ls(drift_time_files[0])

    with PdfPages(args.plot_file) as pdf:
        for det in dets:
            msg = f"... processing {det}"
            log.info(msg)

            # 3. get data observables
            with perf_block("get_data_drift_times()"):
                dts = get_data_drift_times(data, det, ranges=settings.energy_range)

                n = {run: len(dt) for run, dt in dts.items()}

                data_dt_obs, weights, edges = get_data_drift_time_obs(
                    np.concatenate(dts.values()), **settings.dt_kwargs
                )

            with perf_block("plot_drift_time_obs()"):
                fig = plot_drift_time_obs(
                    np.concatenate(dts.values()), data_dt_obs, weights, edges
                )
                decorate(fig)
                pdf.savefig()
                plt.close(fig)

            log.info("... found data observables (%f, %f)", *data_dt_obs)

            # 4. get mc observables
            with perf_block("get_simulated_drift_times()"):
                drift_times_mc, grid_info = get_simulated_drift_times(
                    run_files,
                    det,
                    n,
                    ranges=settings.energy_range,
                )
                log.info("... found MC drift times.")

                depv, slope, dt_obs1, dt_obs2 = get_simulated_drift_time_obs(
                    drift_times_mc, grid_info, **settings.dt_kwargs
                )
                log.info(
                    "... found MC for %d parameters between [%f -- %f] and [%f -- %f]",
                    len(depv),
                    dt_obs1.min(),
                    dt_obs1.max(),
                    dt_obs2.min(),
                    dt_obs2.max(),
                )

            dt_chi2 = get_drift_time_chi2(
                data_dt_obs, (dt_obs1, dt_obs2), settings.drift_time_weight
            )
            log.info(
                "... found chi2 for %d parameters between [%f -- %f]",
                len(depv),
                dt_chi2.min(),
                dt_chi2.max(),
            )
            with perf_block("plot_cost_surface()"):
                fig, _, best_dep, best_slope, best_cost = plot_cost_surface(
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

            out[det] = {
                "slope": float(best_slope),
                "depletion_voltage": float(best_dep),
            }

    dbetto.utils.write_dict(out, args.pars_file)

    print_perf()


if __name__ == "__main__":
    main()
