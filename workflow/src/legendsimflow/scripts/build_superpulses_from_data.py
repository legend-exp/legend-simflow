# ruff: noqa: I002

# Copyright (C) 2026 Giovanna Saleh
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
from pathlib import Path

import numpy as np
from dbetto import AttrsDict
from legendmeta import LegendMetadata
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from legendsimflow import metadata as mutils
from legendsimflow.plot import decorate
from legendsimflow.superpulses import (
    Slice,
    Superpulse,
    compute_chi2,
    get_wfs_for_slice,
    lookup_superpulse_inputs,
    lookup_wfs_indices,
    plot_chi2_cut,
    plot_superpulses,
    plot_wfs_and_superpulse,
    write_superpulses,
)

MIN_NUMBER_WFS = 50
TARGET_WFS = 200
CHI2_THRESHOLD = 3

EVT_TIER_NAME = "pet"

CHARGE_OUTPUT = "wf_pz_win"
CURR_OUTPUT = "curr_av"

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


def main(argv=None):
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

    args = parser.parse_args(argv)

    logger.info(
        "--- Starting Superpulse Generation for %s (%s) ---", args.detector, args.runid
    )

    lmeta = LegendMetadata(str(args.meta))

    raw_files, evt_files, dsp_config, tab_map = lookup_superpulse_inputs(
        l200data=args.l200data,
        metadata=lmeta,
        runid=args.runid,
        hpge=args.detector,
        max_files=args.max_files,
        evt_tier_name=EVT_TIER_NAME,
    )
    file_info = AttrsDict(
        {
            "raw": [str(f) for f in raw_files],
            "evt": [str(f) for f in evt_files],
        }
    )

    logger.info("found %d files per tier", len(file_info.evt))
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

    # extract the indices of waveforms
    # returns a list of indices for each slice
    # for each we have a AttrsDict with fields
    # "hit_idx", "file_idx", "n_sel"

    wf_indices = lookup_wfs_indices(
        slices,
        detector=args.detector,
        evt_files=file_info.evt,
        n_target=TARGET_WFS,
        t0_field="spms/first_t0",
        end_time_field="geds/psd/low_aoe/time",
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    output_lh5 = args.outdir / f"{args.detector}_{period}_{run}_superpulses.lh5"
    output_pdf = args.outdir / f"{args.detector}_{period}_{run}_superpulses.pdf"

    superpulses = {}

    with PdfPages(str(output_pdf)) as pdf:
        for current_slice, slice_wfs_indices in zip(slices, wf_indices, strict=True):
            msg = f"processing {current_slice} ... "
            logger.info(msg)

            if slice_wfs_indices.n_sel < MIN_NUMBER_WFS:
                msg = f"... not enough waveforms {slice_wfs_indices.n_sel} found for {current_slice} skipping"
                logger.warning(msg)
                continue

            # extract the waveform
            wf_data = get_wfs_for_slice(
                file_info.raw,
                f"ch{tab_map[args.detector]}/raw",
                hit_indices=slice_wfs_indices.hit_idx,
                file_indices=slice_wfs_indices.file_idx,
                dsp_config=dsp_config,
                charge_output=CHARGE_OUTPUT,
                current_output=CURR_OUTPUT,
                energy_output="cuspEmax",
                bl_output="bl_std_win",
            )

            if wf_data is None:
                msg = f"... no valid waveforms found for {current_slice} skipping"
                logger.warning(msg)
                continue

            # get a preliminary superpulse
            prelim_sp = Superpulse(
                np.nanmean(wf_data.charge_wfs, axis=0),
                np.nanmean(wf_data.current_wfs, axis=0),
                wf_data.charge_times,
                wf_data.current_times,
                slice=current_slice,
                detector=args.detector,
                n_events_final=wf_data.charge_wfs.shape[0],
                n_events_preliminary=wf_data.charge_wfs.shape[0],
            )

            # chi2 cut
            chi2_values = compute_chi2(
                wf_data.charge_wfs,
                prelim_sp,
                bl_std=wf_data.bl_std,
                cuspEmax=wf_data.energy,
            )

            # cut waveforms above threshold
            sel_charge_wfs = wf_data.charge_wfs[chi2_values < CHI2_THRESHOLD]
            sel_current_wfs = wf_data.current_wfs[chi2_values < CHI2_THRESHOLD]

            n_sel = len(sel_charge_wfs)

            final_sp = Superpulse(
                np.nanmean(sel_charge_wfs, axis=0),
                np.nanmean(sel_current_wfs, axis=0),
                wf_data.charge_times,
                wf_data.current_times,
                slice=current_slice,
                detector=args.detector,
                n_events_final=n_sel,
                n_events_preliminary=wf_data.charge_wfs.shape[0],
            )

            if n_sel < MIN_NUMBER_WFS:
                msg = f"... not enough waveforms {n_sel} selected for {current_slice} skipping"
                logger.warning(msg)
                continue

            superpulses[current_slice] = final_sp

            # Plots
            fig, _ = plot_wfs_and_superpulse(
                wf_data.charge_times,
                wf_data.current_times,
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
                wf_data.charge_times,
                wf_data.charge_wfs,
                final_sp,
                "charge",
            )

            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

            fig, _ = plot_chi2_cut(
                chi2_values,
                CHI2_THRESHOLD,
                wf_data.current_times,
                wf_data.current_wfs,
                final_sp,
                "current",
            )

            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

        # Write superpulses to disk
        if len(superpulses) == 0:
            msg = "No superpulses have been produced!"
            raise RuntimeError(msg)

        msg = f"writing to {output_lh5} ..."
        logger.info(msg)
        write_superpulses(superpulses, str(output_lh5), args.detector, wo_mode="of")

        # Plot superpulses comparison
        fig, _ = plot_superpulses(
            str(output_lh5), args.detector, curve="charge", xlim=(-2000, 1000)
        )
        decorate(fig)
        pdf.savefig(fig)
        plt.close(fig)

        # plot current superpulses
        fig, _ = plot_superpulses(str(output_lh5), args.detector, curve="current")
        decorate(fig)
        pdf.savefig(fig)
        plt.close(fig)

    logger.info("summary plots saved to %s", args.outdir)
    logger.info("done!")

    return 0


if __name__ == "__main__":
    main()
