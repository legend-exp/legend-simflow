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
from matplotlib.backends.backend_pdf import PdfPages
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import utils
from legendsimflow.impurity_tuning import (
    get_drift_time,
    get_drift_time_chi2,
    get_drift_time_obs,
    get_drift_time_obs_mc,
    get_drift_times_mc,
    get_run_mapping,
    plot_drift_time_obs,
    plot_surface,
    read_data,
)
from legendsimflow.metadata import get_simconfig
from legendsimflow.plot import decorate
from legendsimflow.scripts import log_script_invocation

DEFAULT_SETTINGS = {
    "drift_time_weight": 50,  # ns
    "wf_weight": 0.5,  # arb
    "dt_kwargs": {"percentile": 90, "smoothing": 50, "peak_threshold": 0.25},
}


@snakemake_compatible(
    mapping={
        "elecmod": "input.elecmod",
        "drift_time_files": "input.drift_time",
        "data_path": "input.data_path",
        "pars_file": "output.pars_file",
        "plot_file": "output.plot_file",
        "settings": "input.settings",
        "runids": "params.runids",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    """Extract the HPGe electronics model for a LEGEND run.

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
        "--runids",
        required=True,
        nargs="+",
        help="list of runids to process.",
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

    settings = (
        dbetto.AttrsDict(dbetto.utils.load_dict(args.settings))
        if args.settings is not None
        else dbetto.AttrsDict(DEFAULT_SETTINGS)
    )
    log_script_invocation(log, "extract-hpge-impurity-model", parser, args)

    # 1. load data
    msg = f"... loading data from runs {args.runids} and {args.data_path}"
    log.info(msg)
    data = read_data(args.data_path, args.runids)

    simid_mapping = get_run_mapping(
        get_simconfig(config, "hit", simid=None), args.runids
    )

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

            n = {run: len(dt) for run, dt in dts.items()}

            data_dt_obs, weights, edges = get_drift_time_obs(
                np.concatenate(dts.values()), **settings.dt_kwargs
            )

            fig = plot_drift_time_obs(
                np.concatenate(dts.values()), data_dt_obs, weights, edges
            )
            decorate(fig)
            pdf.savefig()
            plt.close(fig)

            log.info("... found data observables (%f, %f)", *data_dt_obs)

            # 4. get mc observables
            weights, dt_mc = get_drift_times_mc(
                args.drift_time_files, det, simid_mapping, n
            )

            depv, slope, dt_obs1, dt_obs2 = get_drift_time_obs_mc(
                dt_mc, grid_info, weights, **settings.dt_kwargs
            )
            dt_chi2 = get_drift_time_chi2(
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

            out[det] = {
                "slope": float(best_slope),
                "depletion_voltage": float(best_dep),
            }

    dbetto.utils.write_dict(out, args.pars_file)


if __name__ == "__main__":
    main()
