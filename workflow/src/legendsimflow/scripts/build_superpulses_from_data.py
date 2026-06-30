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
from pathlib import Path

import dbetto
import hist
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
import numpy as np
from dbetto import AttrsDict
from legendmeta import LegendMetadata
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import utils
from legendsimflow.metadata import get_par_settings
from legendsimflow.plot import decorate
from legendsimflow.scripts import log_script_invocation
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

DEFAULT_SETTINGS = {
    "min_number_wfs": 10,
    "target_wfs": 100,
    "chi2_threshold": 3,
    "max_files": None,
    "evt_tier_name": "pet",
    "charge_output": "wf_pz_win",
    "curr_output": "curr_av",
    "energy_output": "cuspEmax",
    "drift_time_slices": "1000:200:2000",
    "t0_field": "spms/event_t0",
    "end_time_field": "geds/psd/low_aoe/time",
}


def plot_drift_times(dts, slices, tit):
    max_dt = max(cslice.drift_time_range[1] for cslice in slices)

    h = hist.new.Reg(1000, 0, max_dt + 100).Double().fill(dts)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.set_xlabel("Drift time [ns]")
    ax.set_ylabel("counts")

    h.plot(ax=ax, yerr=False, histtype="step")

    for cslice in slices:
        ax.axvline(cslice.drift_time_range[0], linestyle="--", color="grey")
    ax.axvline(slices[-1].drift_time_range[1], linestyle="--", color="grey")
    ax.set_title(tit)
    return fig, ax


@snakemake_compatible(
    mapping={
        "detector": "wildcards.hpge_detector",
        "runid": "params.runids",
        "meta": "config.paths.metadata",
        "l200data": "config.paths.l200data",
        "output_file": "output.superpulses",
        "plot_file": "output.plots",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build HPGe data superpulses from LEGEND-200 data."
    )
    parser.add_argument(
        "--detector", type=str, required=True, help="Detector name (e.g., V03422A)"
    )
    parser.add_argument(
        "--runid",
        type=str,
        required=True,
        help="One or more run IDs (e.g., l200-p16-r006-ssc)",
    )
    parser.add_argument(
        "--meta", type=Path, default=None, help="Path to legend-metadata"
    )
    parser.add_argument(
        "--l200data",
        type=Path,
        default=None,
        help="Path to L200 data production directory",
    )
    parser.add_argument(
        "--output-file", required=True, help="Path to save output files"
    )
    parser.add_argument("--plot-file", required=True, help="Path to save plot files.")

    parser.add_argument("--log-file", default=None, help="log file")

    parser.add_argument(
        "--simflow-config",
        "--config",
        dest="simflow_config",
        required=True,
        help="simflow config YAML path",
    )
    args = parser.parse_args()

    plot_file = Path(args.plot_file)
    output_lh5 = Path(args.output_file)
    plot_file.parent.mkdir(parents=True, exist_ok=True)
    output_lh5.parent.mkdir(parents=True, exist_ok=True)

    simflow_config = utils.init_simflow_context(
        args.simflow_config, workflow=None
    ).config
    log = ldfs.utils.build_log(
        simflow_config.metadata.simprod.config.logging, args.log_file
    )
    log_script_invocation(log, "build-superpulses-from-data", parser, args)

    settings = DEFAULT_SETTINGS | get_par_settings(simflow_config, "superpulses")
    settings = dbetto.AttrsDict(settings)

    meta = (
        Path(args.meta)
        if args.meta is not None
        else Path(simflow_config.paths.metadata)
    )

    l200data = (
        Path(args.l200data)
        if args.l200data is not None
        else getattr(simflow_config.paths, "l200data", None)
    )
    if l200data is None:
        msg = "l200data is not configured and no --l200data path was provided"
        raise RuntimeError(msg)

    runids = args.runid
    runids = runids.strip("[]")
    runids = [r.strip(" '\"") for r in runids.split(",")] if "," in runids else [runids]

    log.info(
        "building superpulses for %s from %s runs (%d in total)",
        args.detector,
        runids,
        len(runids),
    )

    lmeta = LegendMetadata(str(meta))
    raw_files = []
    evt_files = []
    dsp_config = None
    tab_map = {}
    seen_data_runids = set()

    for runid in runids:
        raw_run_files, evt_run_files, dsp_cfg_file, tab_map_run, data_runid = (
            lookup_superpulse_inputs(
                l200data=l200data,
                metadata=lmeta,
                runid=runid,
                hpge=args.detector,
                max_files=settings.max_files,
                evt_tier_name=settings.evt_tier_name,
            )
        )

        # several physics runs can resolve to the same reference calibration
        # run; collect each set of input files only once to avoid duplicating
        # events in the superpulse (and reprocessing the same files)
        if data_runid in seen_data_runids:
            log.debug("skipping %s: data already collected via %s", runid, data_runid)
            continue
        seen_data_runids.add(data_runid)

        raw_files.extend(raw_run_files)
        evt_files.extend(evt_run_files)

        if dsp_config is None:
            dsp_config = dsp_cfg_file
        elif dsp_config != dsp_cfg_file:
            log.warning(
                "multiple DSP configs found across runs (%s, %s); using %s",
                dsp_config,
                dsp_cfg_file,
                dsp_config,
            )
        tab_map.update(tab_map_run)

    if dsp_config is None:
        msg = f"no superpulse input files found for detector {args.detector}"
        raise RuntimeError(msg)

    file_info = AttrsDict(
        {
            "raw": [str(f) for f in raw_files],
            "evt": [str(f) for f in evt_files],
        }
    )

    log.info(
        "found %d raw files and %d evt files",
        len(file_info.raw),
        len(file_info.evt),
    )
    log.info("using DSP config: %s", dsp_config)

    low, step, high = (int(val) for val in settings.drift_time_slices.split(":"))

    slices = [
        Slice(
            energy_range=(1500.0, 2000.0),
            drift_time_range=(float(dt_start), float(dt_start + step)),
        )
        for dt_start in range(low, high, step)
    ]

    # extract the indices of waveforms
    # returns a list of indices for each slice
    # for each we have a AttrsDict with fields
    # "hit_idx", "file_idx", "n_sel"

    wf_indices, drift_times = lookup_wfs_indices(
        slices,
        detector=args.detector,
        evt_files=file_info.evt,
        n_target=settings.target_wfs,
        t0_field=settings.t0_field,
        end_time_field=settings.end_time_field,
    )

    superpulses = {}

    Path(output_lh5).parent.mkdir(parents=True, exist_ok=True)
    Path(plot_file).parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(str(plot_file)) as pdf:
        fig, _ = plot_drift_times(
            drift_times, slices, tit=f"{args.detector} drift time distribution"
        )
        decorate(fig)
        pdf.savefig()
        plt.close(fig)

        for current_slice, slice_wfs_indices in zip(slices, wf_indices, strict=True):
            msg = f"processing {current_slice} ... "
            log.info(msg)

            if slice_wfs_indices.n_sel < settings.min_number_wfs:
                msg = f"... not enough waveforms {slice_wfs_indices.n_sel} found for {current_slice} skipping"
                log.warning(msg)
                continue

            # extract the waveform
            wf_data = get_wfs_for_slice(
                file_info.raw,
                f"ch{tab_map[args.detector]}/raw",
                hit_indices=slice_wfs_indices.hit_idx,
                file_indices=slice_wfs_indices.file_idx,
                dsp_config=dsp_config,
                charge_output=settings.charge_output,
                current_output=settings.curr_output,
                energy_output=settings.energy_output,
                bl_output="bl_std_win",
            )

            if wf_data is None:
                msg = f"... no valid waveforms found for {current_slice} skipping"
                log.warning(msg)
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
            sel_charge_wfs = wf_data.charge_wfs[chi2_values < settings.chi2_threshold]
            sel_current_wfs = wf_data.current_wfs[chi2_values < settings.chi2_threshold]

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

            if n_sel < settings.min_number_wfs:
                msg = f"... not enough waveforms {n_sel} selected for {current_slice} skipping"
                log.warning(msg)
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
                settings.chi2_threshold,
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
                settings.chi2_threshold,
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
        log.info(msg)
        write_superpulses(superpulses, str(output_lh5), args.detector, wo_mode="of")

        # Plot superpulses comparison
        fig, _ = plot_superpulses(str(output_lh5), args.detector, curve="charge")
        decorate(fig)
        pdf.savefig(fig)
        plt.close(fig)

        # plot current superpulses
        fig, _ = plot_superpulses(str(output_lh5), args.detector, curve="current")
        decorate(fig)
        pdf.savefig(fig)
        plt.close(fig)

    log.info("summary plots saved to %s", plot_file)
    log.info("done!")

    return 0


if __name__ == "__main__":
    main()
