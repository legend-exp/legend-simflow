"""Tune the electronics response parameters (sigma, tau) against data superpulses.
 
Reads an ideal pulse shape library and data superpulses from LH5, fits the
Gaussian sigma and exponential tau of the system response kernel by minimising
the mean RMS between simulated and measured current superpulses, and writes
the best-fit parameters to a JSON file.
 
"""
 
from __future__ import annotations
 
import argparse
import json
import logging
 
from lgdo import lh5
from pathlib import Path  

from matplotlib import pyplot as plt  
from matplotlib.backends.backend_pdf import PdfPages  
from legendsimflow.plot import decorate


from legendsimflow.superpulses import read_superpulses
from legendsimflow.electronics_tuning import (
    fit_electronics_parameters,
    plot_best_fit,
    plot_convergence,
    get_ideal_wfs_in_slices,
)
 
logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)
 
 
def main():
    parser = argparse.ArgumentParser(
        description="Fit electronics response parameters (sigma, tau) "
        "against data superpulses."
    )
 
    # Required I/O
    parser.add_argument(
        "--detector", type=str, required=True, help="Detector name (e.g. V03422A)"
    )
    parser.add_argument(
        "--ideal-lib", type=str, required=True, help="Path to ideal PSL LH5 file"
    )
    parser.add_argument(
        "--superpulses",
        type=str,
        required=True,
        help="Path to data superpulses LH5 file",
    )
    parser.add_argument(
        "--output-file", type=str, required=True, help="Path to output JSON file"
    )
 
    # Optional settings
    parser.add_argument(
        "--angle", type=str, default="000", help="Crystal axis angle tag (default: 000)"
    )
    parser.add_argument(
        "--sigma-start",
        type=float,
        default=10.0,
        help="Initial value for sigma in ns (default: 10)",
    )
    parser.add_argument(
        "--tau-start",
        type=float,
        default=70.0,
        help="Initial value for tau in ns (default: 70)",
    )
    parser.add_argument(
        "--sigma-limits",
        type=float,
        nargs=2,
        default=[0.1, 100.0],
        metavar=("LO", "HI"),
        help="Bounds for sigma in ns (default: 0.1 100)",
    )
    parser.add_argument(
        "--tau-limits",
        type=float,
        nargs=2,
        default=[0.1, 500.0],
        metavar=("LO", "HI"),
        help="Bounds for tau in ns (default: 0.1 500)",
    )
    parser.add_argument(
        "--comparison-window",
        type=float,
        nargs=2,
        default=None,
        metavar=("T_MIN", "T_MAX"),
        help="Comparison window in ns relative to the current peak (default: full overlap)",
    )
    parser.add_argument(
        "--max-calls",
        type=int,
        default=5000,
        help="Maximum Minuit function evaluations (default: 5000)",
    )
    parser.add_argument(
        "--plot-dir",
        type=str,
        default=None,
        help="Directory for diagnostic plots (default: no plots)",
    )
    parser.add_argument(
        "--plot-window",
        type=float,
        nargs=2,
        default=None,
        metavar=("T_MIN", "T_MAX"),
        help="X-axis limits for the best-fit overlay in ns (default: same as comparison window)",
    )
 
    args = parser.parse_args()
 
    # Load inputs
    logger.info("reading ideal library from %s ...", args.ideal_lib)
    ideal_lib = lh5.read(args.detector, args.ideal_lib)
 
    logger.info("reading data superpulses from %s ...", args.superpulses)
    data_superpulses = read_superpulses(args.superpulses, args.detector)
    logger.info("loaded %d drift-time slices", len(data_superpulses))
 
    comparison_window = tuple(args.comparison_window) if args.comparison_window else None
 
    # Prepare ideal waveforms
    logger.info("selecting ideal waveforms per slice ...")
    ideal_wfs = get_ideal_wfs_in_slices(ideal_lib, data_superpulses, angle=args.angle)

    # Run fit
    logger.info("starting fit (sigma0=%.1f, tau0=%.1f) ...", args.sigma_start, args.tau_start)
    result = fit_electronics_parameters(
        **ideal_wfs,
        data_superpulses=data_superpulses,
        sigma_start=args.sigma_start,
        tau_start=args.tau_start,
        sigma_limits=tuple(args.sigma_limits),
        tau_limits=tuple(args.tau_limits),
        comparison_window=comparison_window,
        max_calls=args.max_calls,
    )

    # Print summary
    logger.info("")
    logger.info("=" * 50)
    logger.info("  sigma = %.4f ns", result["sigma"])
    logger.info("  tau   = %.4f ns", result["tau"])
    logger.info("  RMS   = %.6f", result["best_rms"])
    logger.info("=" * 50)
 
    # Write output 
    output = {
        "detector": args.detector,
        "angle": args.angle,
        "sigma": result["sigma"],
        "tau": result["tau"],
        "best_rms": result["best_rms"],
    }
 
    with open(args.output_file, "w") as f:
        json.dump(output, f, indent=2)
    logger.info("results written to %s", args.output_file)
 
    # Plots
    if args.plot_dir is not None:

        plot_dir = Path(args.plot_dir)
        plot_dir.mkdir(parents=True, exist_ok=True)

        plot_window = tuple(args.plot_window) if args.plot_window else None

        with PdfPages(str(plot_dir / "tuning_diagnostics.pdf")) as pdf:
            fig, _ = plot_convergence(result)
            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

            fig, _ = plot_best_fit(
                result,
                data_superpulses,
                comparison_window=comparison_window,
                plot_window=plot_window,
            )
            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

        logger.info("saved tuning_diagnostics.pdf")
 
 
if __name__ == "__main__":
    main()