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


import argparse
import logging
from contextlib import nullcontext
from pathlib import Path

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
import lh5
import reboost
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
    """Fit the electronics response over a grid of pulse-shape simulation parameters.

    The ideal (noise-free) waveform library ``--ideal-lib`` holds, for each
    detector, the groups ``psl_scan/<slope>/<depv>``, one per point of a grid
    of impurity-curve slope and depletion voltage, plus an ``info`` group
    with the grid definition (``slope_min``, ``slope_step``, ``dep_min``,
    ``dep_step``): the group names carry the grid indices, not the physical
    values. ``--superpulses`` and ``--settings`` are the data superpulses and
    the fit configuration used for a single-point fit.

    Every grid point is fitted separately, and the parameters are written to
    the YAML file ``--pars-file`` following the same grid layout: entry
    ``[<slope>][<depv>]`` holds the detector name, the ``angle`` the
    superpulses were taken at, the fitted Gaussian width ``sigma`` and
    exponential time constant ``tau`` in ns, and the residual ``rms`` of the
    fit (plus the A/E of data and simulation if plots are produced). The
    ``info`` entry repeats the grid definition and ``best_fit`` is a copy of
    the point with the smallest ``rms``, with the indices it was found at
    added as ``slope`` and ``depv``. Diagnostic plots, one set per grid point,
    go to ``--plot-file``.
    """
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

    if args.plot_file is not None:
        plot_dir = Path(args.plot_file).parent
        plot_dir.mkdir(parents=True, exist_ok=True)

    perf_block, print_perf, _ = reboost.make_profiler()

    output = {}
    with (
        PdfPages(args.plot_file) if args.plot_file is not None else nullcontext() as pdf
    ):
        for slope_group in lh5.ls(args.ideal_lib, f"{args.hpge_detector}/psl_scan/"):
            slope = slope_group.split("/")[-1]
            output[slope] = {}

            log.debug("... reading ideal waveforms from %s ...", slope)

            for depv_group in lh5.ls(
                args.ideal_lib, f"{args.hpge_detector}/psl_scan/{slope}/"
            ):
                depv = depv_group.split("/")[-1]

                with perf_block("read_ideal_wfs()"):
                    ideal_lib = lh5.read(
                        f"{args.hpge_detector}/psl_scan/{slope}/{depv}", args.ideal_lib
                    )

                    # Prepare ideal waveforms
                    ideal_wfs = get_ideal_wfs_all_slices(
                        ideal_lib,
                        data_superpulses,
                        angle=settings.angle,
                        max_num_superpulses=settings.max_num_superpulses,
                    )

                if not ideal_wfs["ideal_wfs_slice"]:
                    log.warning(
                        "no ideal waveforms matched any data superpulse slice for "
                        "slope %s, depv %s, skipping this scan point",
                        slope,
                        depv,
                    )
                    continue

                # Run fit
                log.info(
                    "starting fit (sigma0=%.1f, tau0=%.1f) ...",
                    settings.sigma_start,
                    settings.tau_start,
                )
                with perf_block("fit_electronics_parameters()"):
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

                # plots
                if pdf is not None:
                    with perf_block("plots()"):
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
                        output[slope][depv]["aoe_data"] = data_amax
                        output[slope][depv]["aoe_mc"] = mc_amax

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

    # get the global best fit pars

    best_rms = float("inf")
    best_pars = None

    for slope, slope_dict in output.items():
        for depv, info in slope_dict.items():
            if info["rms"] < best_rms:
                best_rms = info["rms"]
                best_pars = dict(info)
                best_pars["slope"] = slope
                best_pars["depv"] = depv

    step_info = {
        k: float(v.view_as())
        for k, v in lh5.read(f"{args.hpge_detector}/info", args.ideal_lib).items()
    }

    output["info"] = step_info

    if best_pars is not None:
        output["best_fit"] = best_pars
    else:
        msg = "Something went badly wrong, no best fit parameters found!"
        raise RuntimeError(msg)

    print_perf()

    dbetto.utils.write_dict(output, pars_file)
    log.info("... results written to %s", args.pars_file)


if __name__ == "__main__":
    main()
