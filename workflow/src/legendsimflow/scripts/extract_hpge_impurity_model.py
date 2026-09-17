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

@snakemake_compatible(
    mapping={
        "hpge_detector": "wildcards.hpge_detector",
        "elecmod": "input.elecmod",
        "drift_time": "input.drift_time",
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
        "--hpge-detector",
        required=True,
        help="HPGe detector name",
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

