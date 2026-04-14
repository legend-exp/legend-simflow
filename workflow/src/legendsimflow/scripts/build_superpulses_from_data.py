"""
Script to generate HPGe superpulses using the legend-simflow API.
"""
from __future__ import annotations

import argparse
import glob
import logging
from pathlib import Path

from legendmeta import LegendMetadata
from legendsimflow.hpge_pars import lookup_currmod_fit_inputs
from legendsimflow.superpulses import Slice, build_superpulse_for_slice, write_superpulses_to_lh5

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser()
    
    # Required arguments
    parser.add_argument("--detector", type=str, required=True, 
                        help="Detector name (e.g., V03422A)")
    parser.add_argument("--runid", type=str, required=True, 
                        help="Run ID (e.g., l200-p16-r006-ssc)")
    parser.add_argument("--meta", type=Path, required = True, 
                        help="Path to legend-metadata")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="Directory to save the final LH5 files")
    
    args = parser.parse_args()

    logger.info("--- Starting Superpulse Generation for %s (%s) ---", args.detector, args.runid)

    # ----------------------------------------------------------------------------------------------------------------------------------------------
    # Extract period and run - Example: "l200-p16-r008-ssc" -> period="p16", run="r008"
    try:
        parts = args.runid.split('-')
        period = parts[1]
        run = parts[2]
        run_type = parts[3] # This will be 'ssc', 'phy', or 'cal'
    except IndexError:
        raise ValueError(f"runid '{args.runid}' does not match expected format l200-period-run-type")

    # Data Files
    data_path = Path("/global/cfs/cdirs/m2676/data/lngs/l200/public/prodenv/prod-blind/auto/latest")
    raw_files = sorted(glob.glob(f"{data_path}/../../ref/v3.0.0/generated/tier/raw/ssc/{period}/{run}/*.lh5"))[0:10]
    dsp_files = sorted(glob.glob(f"{data_path}/generated/tier/dsp/ssc/{period}/{run}/*.lh5"))[0:10]
    evt_files = sorted(glob.glob(f"{data_path}/generated/tier/evt/ssc/{period}/{run}/*.lh5"))[0:10]
    
    if not evt_files:
        raise FileNotFoundError(f"No EVT files found matching ssc/{period}/{run}/*.lh5")
    
    logger.info("Found %d files for each tier (limited to first 10)", len(evt_files))

    dsp_fields = ["tp_aoe_max"] # Required DSP fields 

    # Setup Metadata and Channel Map
    lmeta = LegendMetadata(str(args.meta))
    start_key = lmeta.datasets.runinfo[period][run][run_type]["start_key"] 
    chmap = lmeta.channelmap(on = start_key)
    tab_map = {args.detector: chmap[args.detector]["daq"]["rawid"]}

    
    # Retrieve the DSP config
    logger.info("Retrieving DSP configuration...")
    _, dsp_config, _, _ = lookup_currmod_fit_inputs(
        l200data=data_path,
        metadata=lmeta,
        runid=args.runid,
        hpge=args.detector
    )
    logger.info("Found DSP config: %s", dsp_config)

    # Slices (Pixels) to process
    # Fixed energy: 1500-2000 keV. Drift time: 500-2500 ns in 200 ns steps.
    slices = []
    for dt_start in range(500, 2500, 200):
        slices.append(Slice(
            energy_range=(1500.0, 2000.0), 
            drift_time_range=(float(dt_start), float(dt_start + 200))
        ))

    # Build superpulse for each slice
    superpulses = {}
    for current_slice in slices:
        logger.info("\nProcessing %s...", current_slice)
        
        sp = build_superpulse_for_slice(
            slice=current_slice,
            detector=args.detector,
            tab_map=tab_map,
            evt_files=evt_files,
            dsp_files=dsp_files,
            raw_files=raw_files,
            dsp_config=dsp_config,
            dsp_fields=dsp_fields,
            n_target_wfs=100,
            chi2_threshold=3.0
        )
        
        superpulses[current_slice] = sp
        logger.info("Golden waveforms: %d/%d", sp.n_events_final, sp.n_events_preliminary)

    # Write superpulses to disk
    args.outdir.mkdir(parents=True, exist_ok=True)
    output_file = args.outdir / f"{args.detector}_{period}_{run}_superpulses.lh5"
    
    logger.info("\nWriting all slices to %s...", output_file)
    write_superpulses_to_lh5(superpulses, str(output_file), args.detector)
    logger.info("Done!")

if __name__ == "__main__":
    main()