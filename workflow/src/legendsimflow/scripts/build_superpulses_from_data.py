"""Script to generate HPGe superpulses using the legend-simflow API."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from legendmeta import LegendMetadata
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from legendsimflow import metadata as mutils
from legendsimflow.plot import decorate
from legendsimflow.superpulses import (
    Slice,
    accumulate_wfs_for_slice,
    apply_chi2_cut,
    compute_superpulse,
    lookup_superpulse_inputs,
    plot_chi2_cut,
    plot_superpulses,
    plot_wfs_and_superpulse,
    write_superpulses_to_lh5,
)

MIN_NUMBER_WFS = 50
TARGET_WFS = 200
CHI2_THRESHOLD = 3.0

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--detector", type=str, required=True, help="Detector name (e.g., V03422A)"
    )
    parser.add_argument(
        "--runid", type=str, required=True, help="Run ID (e.g., l200-p16-r006-ssc)"
    )
    parser.add_argument(
        "--meta", type=Path, required=True, help="Path to legend-metadata"
    )
    parser.add_argument(
        "--l200data",
        type=Path,
        required=True,
        help="Path to L200 data production directory",
    )
    parser.add_argument(
        "--outdir", type=Path, required=True, help="Directory to save output files"
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=None,
        help="limit number of files per tier (default: all)",
    )

    args = parser.parse_args()

    logger.info(
        "--- Starting Superpulse Generation for %s (%s) ---", args.detector, args.runid
    )

    lmeta = LegendMetadata(str(args.meta))

    raw_files, dsp_files, evt_files, dsp_config, tab_map = lookup_superpulse_inputs(
        l200data=args.l200data,
        metadata=lmeta,
        runid=args.runid,
        hpge=args.detector,
        max_files=args.max_files,
    )

    logger.info("found %d files per tier", len(evt_files))
    logger.info("DSP config: %s", dsp_config)

    _, period_int, run_int, _ = mutils.parse_runid(args.runid)
    period = f"p{period_int:02d}"
    run = f"r{run_int:03d}"

    slices = [
        Slice(
            energy_range=(1500.0, 2000.0),
            drift_time_range=(float(dt_start), float(dt_start + 50)),
        )
        for dt_start in range(900, 1900, 50)
    ]

    # Common arguments for accumulate_wfs_for_slice
    file_args = {
        "detector": args.detector,
        "tab_map": tab_map,
        "evt_files": [str(f) for f in evt_files],
        "dsp_files": [str(f) for f in dsp_files],
        "raw_files": [str(f) for f in raw_files],
        "dsp_config": dsp_config,
        "dsp_fields": ["tp_aoe_max"],
    }

    # Output paths
    args.outdir.mkdir(parents=True, exist_ok=True)
    output_lh5 = args.outdir / f"{args.detector}_{period}_{run}_superpulses.lh5"
    output_pdf = args.outdir / f"{args.detector}_{period}_{run}_superpulses.pdf"

    superpulses = {}

    with PdfPages(str(output_pdf)) as pdf:
        for current_slice in slices:
            msg = f"processing {current_slice} ... "
            logger.info(msg)

            # extract the waveforms
            charge_times, current_times, charge_wfs, current_wfs = (
                accumulate_wfs_for_slice(
                    slice=current_slice, n_target=TARGET_WFS, **file_args
                )
            )

            # get a preliminary superpulse
            prelim_sp = compute_superpulse(
                charge_times,
                current_times,
                charge_wfs,
                current_wfs,
                slice=current_slice,
                detector=args.detector,
                n_events=charge_wfs.shape[0],
            )

            # chi2 cut
            chi2_values, sel_charge_wfs, sel_current_wfs = apply_chi2_cut(
                prelim_sp, charge_wfs, current_wfs, chi2_threshold=CHI2_THRESHOLD
            )

            n_sel = len(sel_charge_wfs)

            final_sp = compute_superpulse(
                charge_times,
                current_times,
                sel_charge_wfs,
                sel_current_wfs,
                slice=current_slice,
                detector=args.detector,
                n_events=n_sel,
            )

            if n_sel < MIN_NUMBER_WFS:
                msg = f"... not enough waveforms {n_sel} found for {current_slice} skipping"
                logger.warning(msg)
                continue

            superpulses[current_slice] = final_sp

            # Plots
            fig, _ = plot_wfs_and_superpulse(
                charge_times,
                current_times,
                sel_charge_wfs,
                sel_current_wfs,
                final_sp,
            )
            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

            fig, _ = plot_chi2_cut(
                chi2_values,
                CHI2_THRESHOLD,
                charge_times,
                charge_wfs,
                final_sp,
                "charge",
            )

            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

            fig, _ = plot_chi2_cut(
                chi2_values,
                CHI2_THRESHOLD,
                current_times,
                current_wfs,
                final_sp,
                "current",
            )
            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

        # Write superpulses to disk
        if superpulses:
            msg = f"writing to {output_lh5} ..."
            logger.info(msg)
            write_superpulses_to_lh5(superpulses, str(output_lh5), args.detector)
        else:
            logger.warning("no superpulses produced — nothing to write.")

        # Plot superpulses comparison
        fig, _ = plot_superpulses(
            str(output_lh5), args.detector, curve="charge", xlim=(-2000, 1000)
        )
        pdf.savefig(fig)
        plt.close(fig)

        # plot current superpulses
        fig, _ = plot_superpulses(str(output_lh5), args.detector, curve="current")
        pdf.savefig(fig)
        plt.close(fig)

    logger.info("summary plots saved to %s", args.outdir)
    logger.info("done!")


if __name__ == "__main__":
    main()
