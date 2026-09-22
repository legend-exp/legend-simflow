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

from legendsimflow import nersc, utils
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
        "ideal_psl_scan": "input.ideal_psl_scan",
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

    The ideal (noise-free) waveform library ``--ideal-psl-scan`` holds one
    group per point of a grid of impurity-curve slope and depletion voltage.
    The group names carry the grid indices; the physical values follow from
    the start and the step listed in ``info``::

        <detector>
        |-- psl_scan
        |   |-- slope_0
        |   |   |-- dep_0        # waveform_<angle>_deg, dt, ...
        |   |   `-- dep_1
        |   `-- slope_1
        |       `-- ...
        `-- info                 # slope_min, slope_step, dep_min, dep_step

    ``--superpulses`` and ``--settings`` are the data superpulses and the fit
    configuration of a single-point fit.

    Each grid point is fitted on its own and the results keep the same layout
    in the YAML file ``--pars-file``::

        psl_scan:
            slope_0:
            dep_0:
                detector: V03422A
                angle: "000"       # azimuth of the superpulses, in degrees
                sigma: 12.3        # Gaussian width, in ns
                tau: 47.1          # exponential time constant, in ns
                rms: 0.0021        # residual of the fit
                aoe_data: 1.4      # A/E of data and simulation, only with plots
                aoe_mc: 1.3
            dep_1: ...
        info: ...              # copy of the grid definition above
        best_fit:              # copy of the point with the smallest rms,
          ...                  # with the indices it was found at
          slope: slope_1
          depv: dep_0

    Diagnostic plots, one set per grid point, go to ``--plot-file``.
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
        "--ideal-psl-scan",
        type=str,
        default=None,
        required=False,
        help="path to the LH5 file with the ideal PSL scan grid",
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

        # the LH5 inputs are large, read them through the NERSC read-only mount
        ideal_psl_scan = nersc.dvs_ro(config, args.ideal_psl_scan)
        superpulses = nersc.dvs_ro(config, args.superpulses)
    else:
        logging.basicConfig(
            level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
        )
        log = logging.getLogger(__name__)

        ideal_psl_scan = args.ideal_psl_scan
        superpulses = args.superpulses

    log_script_invocation(log, "extract-hpge-elecmod-scan", parser, args)

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
        superpulses,
    )

    log.info("... reading data superpulses from %s ...", superpulses)
    data_superpulses = read_superpulses(
        superpulses, args.hpge_detector, dt_range_tuning=settings.dt_range_tuning
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

    psl_scan = {}
    with (
        PdfPages(args.plot_file) if args.plot_file is not None else nullcontext() as pdf
    ):
        for slope_group in lh5.ls(ideal_psl_scan, f"{args.hpge_detector}/psl_scan/"):
            slope = slope_group.split("/")[-1]
            psl_scan[slope] = {}

            log.debug("... reading ideal waveforms from %s ...", slope)

            for depv_group in lh5.ls(
                ideal_psl_scan, f"{args.hpge_detector}/psl_scan/{slope}/"
            ):
                depv = depv_group.split("/")[-1]

                with perf_block("read_ideal_wfs()"):
                    ideal_psl = lh5.read(
                        f"{args.hpge_detector}/psl_scan/{slope}/{depv}",
                        ideal_psl_scan,
                    )

                    # Prepare ideal waveforms
                    ideal_wfs = get_ideal_wfs_all_slices(
                        ideal_psl,
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
                psl_scan[slope][depv] = {
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
                        psl_scan[slope][depv]["aoe_data"] = data_amax
                        psl_scan[slope][depv]["aoe_mc"] = mc_amax

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

    for slope, slope_dict in psl_scan.items():
        for depv, info in slope_dict.items():
            if info["rms"] < best_rms:
                best_rms = info["rms"]
                best_pars = dict(info)
                best_pars["slope"] = slope
                best_pars["depv"] = depv

    step_info = {
        k: float(v.view_as())
        for k, v in lh5.read(f"{args.hpge_detector}/grid_info", ideal_psl_scan).items()
    }

    output = {"psl_scan": psl_scan, "grid_info": step_info}

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
