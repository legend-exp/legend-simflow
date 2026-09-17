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
from pathlib import Path

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
import lh5
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import metadata as mutils
from legendsimflow import utils
from legendsimflow.hpge_electronics_tuning import (
    fit_electronics_parameters,
    get_ideal_wfs_all_slices,
    plot_best_fit,
    plot_convergence,
)
from legendsimflow.plot import decorate
from legendsimflow.scripts import log_script_invocation
from legendsimflow.superpulses import (
    plot_current_superpulses_fwhm_and_amplitude,
    read_superpulses,
)

DEFAULT_SETTINGS = {
    "drift_time_weight": 50, # ns
    "wf_weight": 0.5, # arb
}


def get_drift_time_obs_mc(mc):
    depv = []
    obs1 = []
    obs2 = []
    slope = []
    
    for dep in out.fields:
        if dep == "weight":
            continue
        depf = float(dep)
        for s in out[dep].fields:
            slopef = float(s)
            obs = drift_time.drift_time_observables(out[dep][s],weights = out.weight)
    
        
            depv.append(depf)
            slope.append(slopef)    
            obs1.append(obs[0])
            obs2.append(obs[1]-obs[0])

    depv= np.array(depv)
    slope = np.array(slope)
    obs1 = np.array(obs1)
    obs2 = np.array(obs2)
          
    return depv,slope,obs1,obs2

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
        required=True,
        help="input YAML file for the electronics model",
    )
    parser.add_argument(
        "--drift-time",
        required=True,
        help="input LH5 file for the drift times.",
    )
    parser.add_argument(
        "--pars-file",
        required=True,
        help="output YAML file for the electronics model parameters",
    )
    parser.add_argument(
        "--data-files",
        required=True,
        nargs="+",
        help="input LH5 files for the data.",
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

    log = ldfs.utils.build_log(metadata.simprod.config.logging, args.log_file)
    log_script_invocation(log, "extract-hpge-elecmod", parser, args)


    # 1. load data
    data = load_data()

    # 2. get run norms (e.g. from livetime)
    run_norms = get_run_norms(args.data_files, metadata)

    out = {}
    with PdfPages(args.plot_file) as pdf:
        for det in dets:
            
            # 3. get data observables
            data_dt_obs = get_drift_time_obs(data,**settings.dt_kwargs)

            # 4. get mc observables
            
            dt_obs_mc = get_drift_time_obs_mc(drift_times,**settings.dt_kwargs)

            dt_chi2 = get_dt_chi2(data_dt_obs, dt_obs_mc,settings.drift_time_weight)

            # 5. extract electronics model parameters

            elecmod = dbetto.utils.load_dict(args.elecmod)

            wf_chi2 = get_wf_chi2(elecmod,settings.wf_weight)

            # find the best fit
            best_slope, best_dep  = fit_impurities(det,wf_chi2,dt_chi2,pdf)

            out[det] = {"slope": best_slope, "depletion_voltage": best_dep}

    dbetto.utils.write_dict(out, args.pars_file)

if __name__ == "__main__":
    main()



