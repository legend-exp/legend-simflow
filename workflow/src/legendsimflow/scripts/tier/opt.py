# ruff: noqa: I002

# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>,
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

import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
import pyg4ometry
import pygeomtools
import reboost
import reboost.spms
from dbetto import AttrsDict
from dbetto.utils import load_dict
from lh5 import LH5Iterator
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, optics, utils
from legendsimflow.exceptions import SimflowConfigError
from legendsimflow.metadata import get_tier_settings
from legendsimflow.scripts import log_script_invocation


@snakemake_compatible(
    mapping={
        "stp_file": "input.stp_file",
        "optmap_lar": lambda snakemake: (
            snakemake.input.optmap_lar[0] if snakemake.input.optmap_lar else None
        ),
        "geom_file": "input.geom",
        "simstat_part_file": "input.simstat_part_file",
        "usability_file": "input.usability",
        "jobid": "wildcards.jobid",
        "opt_file": "output[0]",
        "log_file": "log[0]",
        "optmap_per_sipm": "params.optmap_per_sipm",
        "scintillator_volume_name": "params.scintillator_volume_name",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(description="Build the opt tier.")
    parser.add_argument("--stp-file", required=True, help="input stp tier file")
    parser.add_argument(
        "--optmap-lar",
        default=None,
        help="LAr optical map file, needed with light_source: optmap",
    )
    parser.add_argument("--geom-file", required=True, help="GDML geometry file")
    parser.add_argument(
        "--simstat-part-file",
        required=True,
        help="simulation statistics partition file",
    )
    parser.add_argument(
        "--usability-file",
        required=True,
        help="detector usability YAML file",
    )
    parser.add_argument("--jobid", required=True, help="job ID wildcard")
    parser.add_argument("--opt-file", required=True, help="output opt tier file")
    parser.add_argument("--log-file", default=None, help="log file")
    parser.add_argument(
        "--optmap-per-sipm",
        action="store_true",
        default=False,
        help="sample photoelectrons per SiPM (default: all SiPMs together)",
    )
    parser.add_argument(
        "--scintillator-volume-name",
        required=True,
        help="name of the scintillator sensitive volume in the geometry",
    )
    parser.add_argument(
        "--simflow-config",
        "--config",
        dest="simflow_config",
        required=True,
        help="simflow config YAML path",
    )
    args = parser.parse_args()

    config = utils.init_simflow_context(args.simflow_config, workflow=None).config

    stp_file = nersc.dvs_ro(config, args.stp_file)
    jobid = args.jobid
    opt_file = args.opt_file
    gdml_file = nersc.dvs_ro(config, args.geom_file)
    log_file = args.log_file
    metadata = config.metadata
    optmap_per_sipm = args.optmap_per_sipm
    scintillator_volume_name = args.scintillator_volume_name
    simstat_part_file = nersc.dvs_ro(config, args.simstat_part_file)
    usability_map = AttrsDict(load_dict(nersc.dvs_ro(config, args.usability_file)))

    opt_file, move2cfs = nersc.make_on_scratch(config, opt_file)

    tier_opt_settings = get_tier_settings(config, "opt")
    buffer_len = tier_opt_settings.buffer_len
    light_source = tier_opt_settings.light_source

    settings_block = f"simprod.config.tier.opt.{config.experiment}.settings"
    if light_source not in ("optmap", "tracked_photons"):
        msg = (
            f"light_source must be 'optmap' or 'tracked_photons', not {light_source!r}"
        )
        raise SimflowConfigError(msg, settings_block)
    if light_source == "optmap":
        if args.optmap_lar is None:
            msg = "--optmap-lar is required with light_source: optmap"
            raise SimflowConfigError(msg, settings_block)
        optmap_lar = nersc.dvs_ro(config, args.optmap_lar)
    else:
        optmap_lar = None

    # setup logging
    log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "tier-opt", parser, args)
    perf_block, print_perf, print_perf_last = reboost.make_profiler()

    # load the geometry and retrieve registered sensitive volume tables
    geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
    scintillators = pygeomtools.detectors.get_all_senstables(
        geom, type_filter="scintillator"
    )
    sipms = pygeomtools.detectors.get_all_senstables(geom, type_filter="optical")

    # fail early and loudly: without a matching volume the loop below would
    # silently do nothing and the output file would never be created. tracked
    # photons do not need the scintillator table
    if light_source == "optmap" and scintillator_volume_name not in scintillators:
        msg = (
            f"scintillator volume {scintillator_volume_name} not found in "
            f"{gdml_file}. scintillator volumes in the geometry: "
            f"{sorted(scintillators)}"
        )
        raise SimflowConfigError(msg, settings_block)

    partitions = load_dict(simstat_part_file)[f"job_{jobid}"]

    # load TCM, to be used to chunk the event statistics according to the run partitioning
    msg = "loading TCM"
    log.debug(msg)
    tcm = lh5.read_as("tcm", stp_file, library="ak")

    # in combined mode there is a single map, so pre-load it once for a little
    # speed up. in per-SiPM mode the maps are loaded on demand, one at a time
    if light_source == "optmap" and not optmap_per_sipm:
        optmap_lar = reboost.spms.load_optmap(optmap_lar, "all")

    # loop over the partitions for this file
    for runid_idx, (runid, evt_idx_range) in enumerate(partitions.items()):
        msg = (
            f"processing partition corresponding to {runid} "
            f"[{runid_idx + 1}/{len(partitions)}], event range {evt_idx_range}"
        )
        log.info(msg)

        if light_source == "tracked_photons":
            # as in the hit tier, each row of a SiPM table gives one output row
            for sipm in sorted(sipms):
                sipm_uid = sipms[sipm].uid
                i_start, n_entries = reboost.get_rows_in_event_range(
                    tcm, sipm_uid, *evt_idx_range
                )
                if n_entries == 0:
                    continue

                usability = optics.get_usability(usability_map, runid, sipm)

                msg = f"processing the tracked photons of SiPM {sipm}"
                log.debug(msg)

                photons = LH5Iterator(
                    stp_file,
                    f"stp/{sipm}",
                    i_start=i_start,
                    n_entries=n_entries,
                    buffer_len=buffer_len,
                    field_mask=["evtid", "t0", "time"],
                )
                optics.write_sipm_photoelectrons(
                    photons,
                    opt_file,
                    sipm,
                    sipm_uid,
                    runid,
                    usability,
                    settings=tier_opt_settings,
                    perf_block=perf_block,
                )
                print_perf_last()
            continue

        # convolve the scintillator table with the optical maps
        scint_uid = scintillators[scintillator_volume_name].uid
        msg = (
            f"looking for data from sensitive volume {scintillator_volume_name} "
            f"table (uid={scint_uid})..."
        )
        log.debug(msg)

        if f"stp/{scintillator_volume_name}" not in lh5.ls(stp_file, "stp/*"):
            msg = (
                f"detector {scintillator_volume_name} not found in {stp_file}. "
                "possibly because it was not read-out or there were no hits recorded"
            )
            log.warning(msg)
            continue

        log.info("processing the 'lar' scintillator table...")

        msg = "looking for indices of hit table rows to read..."
        log.debug(msg)
        i_start, n_entries = reboost.get_rows_in_event_range(
            tcm, scint_uid, *evt_idx_range
        )

        # output table, its UID and the usability of the SiPM
        if optmap_per_sipm:
            targets = [
                (
                    sipm,
                    sipms[sipm].uid,
                    optics.get_usability(usability_map, runid, sipm),
                )
                for sipm in sorted(sipms)
            ]
        else:
            targets = [("all", scint_uid, "on")]

        for sipm, uid, usability in targets:
            msg = f"applying optical map for SiPM {sipm}"
            log.debug(msg)

            steps = LH5Iterator(
                stp_file,
                f"stp/{scintillator_volume_name}",
                i_start=i_start,
                n_entries=n_entries,
                buffer_len=buffer_len,
            )
            stats = optics.write_sipm_photoelectrons(
                steps,
                opt_file,
                sipm,
                uid,
                runid,
                usability,
                settings=tier_opt_settings,
                optmap=optmap_lar,
                perf_block=perf_block,
            )
            optics.log_optmap_coverage(stats, sipm, scintillator_volume_name)
            print_perf_last()

    # a SiPM without photons still gets a table: the evt tier takes the list of
    # channels from the opt file, which must not change from job to job
    if light_source == "tracked_photons":
        tables = lh5.ls(opt_file, "hit/*") if Path(opt_file).exists() else []
        for sipm in sorted(sipms):
            if f"hit/{sipm}" in tables:
                continue
            optics.write_sipm_photoelectrons(
                [optics.empty_photon_table()],
                opt_file,
                sipm,
                sipms[sipm].uid,
                runid,
                "on",
                settings=tier_opt_settings,
            )

    log.debug("building the TCM")
    reboost.build_remage_tcm(opt_file, opt_file)

    with perf_block("move_to_cfs()"):
        move2cfs()

    print_perf()


if __name__ == "__main__":
    main()
