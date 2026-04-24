"""Script to generate HPGe superpulses using the legend-simflow API."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from legendmeta import LegendMetadata

from legendsimflow import metadata as mutils
from legendsimflow import utils
from legendsimflow.superpulses import (
    Slice,
    build_superpulse_for_slice,
    write_superpulses_to_lh5,
)

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


def lookup_superpulse_inputs(  ### I will move this to the superpulses module or to a utility module if the content is ok
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
        "--outdir", type=Path, required=True, help="Directory to save output LH5 files"
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
            drift_time_range=(float(dt_start), float(dt_start + 200)),
        )
        for dt_start in range(1100, 2300, 200)
    ]

    superpulses = {}
    for current_slice in slices:
        logger.info("processing %s...", current_slice)
        sp = build_superpulse_for_slice(
            slice=current_slice,
            detector=args.detector,
            tab_map=tab_map,
            evt_files=[str(f) for f in evt_files],
            dsp_files=[str(f) for f in dsp_files],
            raw_files=[str(f) for f in raw_files],
            dsp_config=dsp_config,
            dsp_fields=["tp_aoe_max"],
            n_target_wfs=100,
            chi2_threshold=3.0,
        )
        superpulses[current_slice] = sp
        logger.info(
            "golden waveforms: %d/%d", sp.n_events_final, sp.n_events_preliminary
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    output_file = args.outdir / f"{args.detector}_{period}_{run}_superpulses.lh5"
    logger.info("writing to %s...", output_file)
    write_superpulses_to_lh5(superpulses, str(output_file), args.detector)
    logger.info("done!")


if __name__ == "__main__":
    main()
