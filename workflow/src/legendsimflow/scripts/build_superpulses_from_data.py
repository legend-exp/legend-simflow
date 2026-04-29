"""Script to generate HPGe superpulses using the legend-simflow API."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
from legendmeta import LegendMetadata
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from legendsimflow import metadata as mutils
from legendsimflow import utils
from legendsimflow.superpulses import (
    Slice,
    accumulate_wfs_for_slice,
    apply_chi2_cut,
    compute_chi2_vs_superpulse,
    compute_superpulse,
    plot_chi2_cut,
    plot_superpulses,
    plot_wfs_and_superpulse,
    trim_and_stack,
    write_superpulses_to_lh5,
)

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


def lookup_superpulse_inputs(
    l200data: str | Path,
    metadata: LegendMetadata,
    runid: str,
    hpge: str,
    max_files: int | None = None,
) -> tuple[list[Path], list[Path], list[Path], Path, dict[str, int]]:
    """Look up all inputs needed to build superpulses for one detector and run.

    Parameters
    ----------
    l200data
        Path to the L200 data production directory.
    metadata
        The metadata instance.
    runid
        LEGEND run identifier, e.g. ``"l200-p16-r008-ssc"``.
    hpge
        Detector name, e.g. ``"V03422A"``.
    max_files
        Limit the number of files per tier. Default: all.

    Returns
    -------
    raw_files, dsp_files, evt_files
        Sorted lists of LH5 file paths for each tier.
    dsp_cfg_file
        Path to the DSP configuration file.
    tab_map
        Mapping ``{detector_name: rawid}``.
    """
    if isinstance(l200data, str):
        l200data = Path(l200data)

    _, period_int, run_int, data_type = mutils.parse_runid(runid)
    period = f"p{period_int:02d}"
    run = f"r{run_int:03d}"
    df_cfg = utils.lookup_dataflow_config(l200data).paths

    raw_files = sorted((df_cfg.tier_raw / data_type / period / run).glob("*.lh5"))[
        :max_files
    ]
    dsp_files = sorted((df_cfg.tier_dsp / data_type / period / run).glob("*.lh5"))[
        :max_files
    ]
    evt_files = sorted((df_cfg.tier_evt / data_type / period / run).glob("*.lh5"))[
        :max_files
    ]

    if not evt_files:
        msg = f"no EVT files found for {runid}"
        raise FileNotFoundError(msg)

    dsp_cfg_regex = r"l200-*-r%-T%-ICPC-dsp_proc_chain.*"
    dsp_cfg_files = list(
        (l200data / "inputs/dataprod/config/tier_dsp").glob(dsp_cfg_regex)
    )
    if dsp_cfg_files == []:
        dsp_cfg_files = list(
            (l200data / "inputs/dataprod/config/tier/dsp").glob(dsp_cfg_regex)
        )
    if len(dsp_cfg_files) != 1:
        msg = f"could not find a suitable dsp config file in {l200data} (or multiple found)"
        raise RuntimeError(msg)

    tstamp = mutils.runinfo(metadata, runid).start_key
    chmap = metadata.channelmap(tstamp)
    tab_map = {hpge: chmap[hpge]["daq"]["rawid"]}

    return raw_files, dsp_files, evt_files, dsp_cfg_files[0], tab_map


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

    # Configuration
    n_target_wfs = 100
    chi2_threshold = 3.0
    _ACCUMULATION_FACTOR = 1.5

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
            logger.info("processing %s...", current_slice)
            file_idx, final_sp = 0, None
            b_c_times, b_i_times, b_c_wfs, b_i_wfs, all_bl, all_ce = (
                [],
                [],
                [],
                [],
                [],
                [],
            )

            while file_idx < len(evt_files):
                n_need = max(
                    1,
                    int(_ACCUMULATION_FACTOR * n_target_wfs)
                    - sum(len(c) for c in b_c_wfs),
                )

                ct, it, cw, iw, bl, ce, file_idx = accumulate_wfs_for_slice(
                    slice=current_slice,
                    start_file_idx=file_idx,
                    n_target=n_need,
                    **file_args,
                )

                if cw is not None:
                    b_c_times.append(ct)
                    b_i_times.append(it)
                    b_c_wfs.append(cw)
                    b_i_wfs.append(iw)
                    if bl is not None:
                        all_bl.append(bl)
                    if ce is not None:
                        all_ce.append(ce)

                if not b_c_wfs:
                    continue

                # Re-trim all batches and stack using the helper
                charge_times, charge_stack = trim_and_stack(b_c_times, b_c_wfs)
                current_times, current_stack = trim_and_stack(b_i_times, b_i_wfs)

                bl_std_stack = np.concatenate(all_bl) if all_bl else None
                cuspEmax_stack = np.concatenate(all_ce) if all_ce else None
                n_preliminary = charge_stack.shape[0]

                # preliminary superpulse
                prelim_sp = compute_superpulse(
                    charge_times,
                    current_times,
                    charge_stack,
                    current_stack,
                    current_slice,
                    args.detector,
                    n_preliminary,
                )

                # chi2
                chi2_values = compute_chi2_vs_superpulse(
                    charge_stack,
                    prelim_sp,
                    bl_std=bl_std_stack,
                    cuspEmax=cuspEmax_stack,
                )

                golden_charge, golden_current, _ = apply_chi2_cut(
                    charge_stack,
                    current_stack,
                    chi2_values,
                    threshold=chi2_threshold,
                )
                n_golden = golden_charge.shape[0]

                logger.info(
                    "golden waveforms: %d/%d (target: %d)",
                    n_golden,
                    n_preliminary,
                    n_target_wfs,
                )

                if n_golden >= n_target_wfs or file_idx >= len(evt_files):
                    if n_golden == 0:
                        logger.warning(
                            "no waveforms survived chi2 cut; "
                            "using preliminary superpulse"
                        )
                        final_sp = prelim_sp
                        golden_charge = charge_stack
                        golden_current = current_stack
                        break

                    if n_golden < n_target_wfs:
                        logger.warning(
                            "all files exhausted: %d golden (target: %d)",
                            n_golden,
                            n_target_wfs,
                        )

                    # final superpulse
                    final_sp = compute_superpulse(
                        charge_times,
                        current_times,
                        golden_charge,
                        golden_current,
                        current_slice,
                        args.detector,
                        n_preliminary,
                    )
                    break

                logger.info("need more waveforms, continuing...")

            if final_sp is None:
                logger.warning("no waveforms found for %s, skipping", current_slice)
                continue

            superpulses[current_slice] = final_sp

            # Plots
            fig, _ = plot_wfs_and_superpulse(
                charge_times,
                current_times,
                golden_charge,
                golden_current,
                final_sp,
            )
            pdf.savefig(fig)
            plt.close(fig)

            fig, _ = plot_chi2_cut(
                chi2_values,
                chi2_threshold,
                charge_times,
                charge_stack,
                final_sp,
                "charge",
            )
            pdf.savefig(fig)
            plt.close(fig)

            fig, _ = plot_chi2_cut(
                chi2_values,
                chi2_threshold,
                current_times,
                current_stack,
                final_sp,
                "current",
            )

            pdf.savefig(fig)
            plt.close(fig)

        # Write superpulses to disk
        if superpulses:
            logger.info("writing to %s...", output_lh5)
            write_superpulses_to_lh5(superpulses, str(output_lh5), args.detector)
        else:
            logger.warning("no superpulses produced — nothing to write.")

        # Plot superpulses comparison
        fig, _ = plot_superpulses(
            str(output_lh5), args.detector, curve="charge", xlim=(-2000, 1000)
        )
        pdf.savefig(fig)
        plt.close(fig)

        fig, _ = plot_superpulses(str(output_lh5), args.detector, curve="current")
        pdf.savefig(fig)
        plt.close(fig)

    logger.info("summary plots saved to %s", args.outdir)
    logger.info("done!")


if __name__ == "__main__":
    main()
