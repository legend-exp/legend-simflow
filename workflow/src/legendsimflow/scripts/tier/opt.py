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
from legendsimflow.metadata import get_opt_scintillators, get_tier_settings
from legendsimflow.scripts import log_script_invocation


@snakemake_compatible(
    mapping={
        "stp_file": "input.stp_file",
        "optmaps": lambda snakemake: [
            f"{key.removeprefix('optmap_')}={path}"
            for key, path in snakemake.input.items()
            if key.startswith("optmap_")
        ],
        "geom_file": "input.geom",
        "simstat_part_file": "input.simstat_part_file",
        "usability_file": "input.usability",
        "jobid": "wildcards.jobid",
        "opt_file": "output[0]",
        "log_file": "log[0]",
        "optmap_per_sipm": "params.optmap_per_sipm",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(description="Build the opt tier.")
    parser.add_argument("--stp-file", required=True, help="input stp tier file")
    parser.add_argument(
        "--optmap",
        dest="optmaps",
        action="append",
        default=[],
        metavar="SOURCE=FILE",
        help=(
            "optical map of a light source, needed with light_source: optmap "
            "for each source without a light model"
        ),
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
        sources = get_opt_scintillators(config)
        optmaps = {}
        for item in args.optmaps:
            name, _, path = item.partition("=")
            optmaps[name] = nersc.dvs_ro(config, path)
        mapped = sorted(n for n, s in sources.items() if s.model is None)
        if sorted(optmaps) != mapped:
            msg = (
                f"one optical map per light source without a light model is "
                f"needed: sources {mapped}, maps given for {sorted(optmaps)}"
            )
            raise SimflowConfigError(msg, settings_block)
        # a light model takes the place of the map of its source
        optmaps |= {
            name: optics.load_light_model(s.model)
            for name, s in sources.items()
            if s.model is not None
        }

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

    if light_source == "optmap":
        # fail early and loudly: without a matching table the output file
        # would never be created. tracked photons do not need these tables
        for name, source in sources.items():
            if source.table not in scintillators:
                msg = (
                    f"scintillator table {source.table} of light source {name} not "
                    f"found in {gdml_file}. scintillator tables in the geometry: "
                    f"{sorted(scintillators)}"
                )
                raise SimflowConfigError(msg, settings_block)
        uid2source = {scintillators[s.table].uid: name for name, s in sources.items()}

        stp_tables = lh5.ls(stp_file, "stp/*")
        for name, source in sources.items():
            if f"stp/{source.table}" not in stp_tables:
                msg = (
                    f"table {source.table} of light source {name} not found in "
                    f"{stp_file}. possibly because it was not read-out or there "
                    "were no hits recorded"
                )
                log.warning(msg)
        # chunk size in rows, the smallest of the scintillator tables
        rows_per_chunk = min(
            (
                LH5Iterator(
                    stp_file, f"stp/{s.table}", buffer_len=buffer_len
                ).buffer_len
                for s in sources.values()
                if f"stp/{s.table}" in stp_tables
            ),
            default=1,
        )

    partitions = load_dict(simstat_part_file)[f"job_{jobid}"]

    # load TCM, to be used to chunk the event statistics according to the run partitioning
    msg = "loading TCM"
    log.debug(msg)
    tcm = lh5.read_as("tcm", stp_file, library="ak")

    # in combined mode there is a single map per source, so pre-load them once
    # for a little speed up. in per-SiPM mode the maps are loaded on demand
    if light_source == "optmap" and not optmap_per_sipm:
        optmaps = {
            name: reboost.spms.load_optmap(m, "all")
            if isinstance(m, (str, Path))
            else m
            for name, m in optmaps.items()
        }

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
                optics.write_tracked_photoelectrons(
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

        # each event window (TCM row) with a scintillator hit gives one output
        # row, where the light of all its scintillator hits is pooled
        first, last = evt_idx_range
        windows, n_windows = optics.find_event_windows(
            tcm[first : last + 1], uid2source
        )
        msg = (
            f"found {n_windows} event windows with scintillator hits, "
            f"{ {n: len(w) for n, (_, w) in windows.items()} } hits per light source"
        )
        log.info(msg)
        if n_windows == 0:
            continue

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
            # the smallest UID of the light sources, independent of their
            # order in the configuration
            targets = [("all", min(uid2source), "on")]

        for sipm, uid, usability in targets:
            msg = f"applying optical maps and light models for SiPM {sipm}"
            log.debug(msg)

            stats = optics.write_scintillation_photoelectrons(
                stp_file,
                windows,
                n_windows,
                sources,
                optmaps,
                opt_file,
                sipm,
                uid,
                runid,
                usability,
                settings=tier_opt_settings,
                rows_per_chunk=rows_per_chunk,
                perf_block=perf_block,
            )
            for name, source_stats in stats.items():
                optics.log_optmap_coverage(
                    source_stats, f"{name}/{sipm}", sources[name].table
                )
            print_perf_last()

    # a SiPM without photons still gets a table: the evt tier takes the list of
    # channels from the opt file, which must not change from job to job
    if light_source == "tracked_photons":
        tables = lh5.ls(opt_file, "hit/*") if Path(opt_file).exists() else []
        for sipm in sorted(sipms):
            if f"hit/{sipm}" in tables:
                continue
            optics.write_tracked_photoelectrons(
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
