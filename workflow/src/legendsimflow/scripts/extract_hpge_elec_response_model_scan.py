# ruff: noqa: I002

# Copyright (C) 2026 Toby Dixon <toby.dixon.23@ucl.ac.uk>, Giovanna Saleh <giovanna.saleh@phd.unipd.it>,
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

"""Tune the electronics response parameters (sigma, tau) against data superpulses.

Reads an ideal pulse-shape library and data superpulses from LH5, fits the
Gaussian sigma and exponential tau of the system response kernel by minimising
the mean RMS between simulated and measured current superpulses, and writes
the best-fit parameters to a YAML file.

"""

import argparse
import logging
from pathlib import Path

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
import lh5
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from snakemake_argparse_bridge import snakemake_compatible

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
    read_superpulses,
)

DEFAULT_SETTINGS = {
    "angle": "000",
    "sigma_start": 10.0,
    "tau_start": 50.0,
    "sigma_limits": (0.0, 200.0),
    "tau_limits": (0.0, 200.0),
    "comparison_window": (-500.0, 500.0),
    "plot_window": (-600.0, 600.0),
    "weight_power": 2.0,
    "max_calls": 1000,
    "dt_range_tuning": (600.0, 3000.0),
    "max_num_superpulses": 5,
    "truncate_gauss": True,
}


@snakemake_compatible(
    mapping={
        "hpge_detector": "wildcards.hpge_detector",
        "ideal_lib": "input.ideal_psl",
        "superpulses": "input.superpulses",
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
        dest="hpge_detector",
        help="HPGe detector name",
    )
    parser.add_argument(
        "--ideal-lib",
        type=str,
        default=None,
        required=False,
        help="Path to ideal psl file",
    )
    parser.add_argument(
        "--superpulses",
        type=str,
        required=False,
        default=None,
        help="Path to data superpulses file",
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

    if args.simflow_config is not None:
        config = utils.init_simflow_context(args.simflow_config, workflow=None).config
        metadata = config.metadata

        log_config = metadata.simprod.config.logging
        log = ldfs.utils.build_log(log_config, args.log_file)

        log_script_invocation(log, "extract-hpge-elecmod-scan", parser, args)
    else:
        logging.basicConfig(
            level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
        )
        log = logging.getLogger(__name__)

    hpge = args.hpge_detector
    pars_file = args.pars_file

    settings = (
        dbetto.AttrsDict(dbetto.utils.load_dict(args.settings))
        if args.settings is not None
        else dbetto.AttrsDict(DEFAULT_SETTINGS)
    )

    log.info(
        "extracting electronics model from superpulses %s in %s ...",
        hpge,
        args.superpulses,
    )

    log.info("... reading data superpulses from %s ...", args.superpulses)
    data_superpulses = read_superpulses(
        args.superpulses, args.hpge_detector, dt_range_tuning=settings.dt_range_tuning
    )

    if not data_superpulses:
        msg = f"no superpulses found in drift time range [{settings.dt_range_tuning[0]:.0f}, {settings.dt_range_tuning[1]:.0f}] ns"
        raise RuntimeError(msg)

    # loop over slope and depv
    comparison_window = tuple(settings.comparison_window)
    plot_window = tuple(settings.plot_window)

    output = {}
    for slope_group in lh5.ls(args.ideal_lib, f"{args.hpge_detector}/"):
        slope = slope_group.split("/")[-1]
        output[slope] = {}

        log.info("... reading ideal waveforms from %s ...", slope)

        for depv_group in lh5.ls(args.ideal_lib, f"{args.hpge_detector}/{slope}/"):
            depv = depv_group.split("/")[-1]
            log.info("... reading ideal waveforms from %s ...", depv)

            ideal_lib = lh5.read(f"{args.hpge_detector}/{slope}/{depv}", args.ideal_lib)

            # Prepare ideal waveforms
            log.info("... selecting ideal waveforms per slice ...")
            ideal_wfs = get_ideal_wfs_all_slices(
                ideal_lib,
                data_superpulses,
                angle=settings.angle,
                max_num_superpulses=settings.max_num_superpulses,
            )

            if not ideal_wfs["ideal_wfs_slice"].keys():
                msg = "no ideal waveforms matched any data superpulse slice"
                raise RuntimeError(msg)
         

            # Run fit
            log.info(
                "starting fit (sigma0=%.1f, tau0=%.1f) ...",
                settings.sigma_start,
                settings.tau_start,
            )
            result = fit_electronics_parameters(
                **ideal_wfs,
                data_superpulses=data_superpulses,
                sigma_start=settings.sigma_start,
                tau_start=settings.tau_start,
                sigma_limits=tuple(settings.sigma_limits),
                tau_limits=tuple(settings.tau_limits),
                comparison_window=comparison_window,
                weight_power=settings.get("weight_power", 0.0),
                max_calls=settings.max_calls,
            )

            # Write output
            output[slope][depv] = {
                "detector": args.hpge_detector,
                "angle": settings.angle,
                "sigma": result["sigma"],
                "tau": result["tau"],
                "rms": result["best_rms"],
            }

    print(output)
    return 1

    # Plots
    if args.plot_file is not None:
        plot_dir = Path(args.plot_file).parent
        plot_dir.mkdir(parents=True, exist_ok=True)

        with PdfPages(str(args.plot_file)) as pdf:
            fig, _ = plot_convergence(result)
            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

            fig, _, data_amax, mc_amax = plot_best_fit(
                result,
                data_superpulses,
                comparison_window=comparison_window,
                plot_window=plot_window,
                detector_name=args.hpge_detector,
            )
            output["aoe_data"] = data_amax
            output["aoe_mc"] = mc_amax
            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

            fig, _, _, _ = plot_best_fit(
                result,
                data_superpulses,
                comparison_window=comparison_window,
                plot_window=plot_window,
                plot_charge=True,
                detector_name=args.hpge_detector,
            )
            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

        log.info("... saved diagnostic plots to %s", args.plot_file)

    dbetto.utils.write_dict(output, pars_file)
    log.info("... results written to %s", args.pars_file)


if __name__ == "__main__":
    main()
