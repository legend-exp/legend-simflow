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
import re
from collections.abc import Iterable
from pathlib import Path

import awkward as ak
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
import numpy as np
from dbetto.utils import load_dict
from lgdo import Scalar, Struct
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, patterns, utils
from legendsimflow.exceptions import SimflowConfigError
from legendsimflow.metadata import get_tier_settings, parse_runid, runinfo
from legendsimflow.partitioning import primaries_per_run
from legendsimflow.pdf import (
    fill_histograms,
    format_pdf_struct,
    group_detectors,
    merge_pdf_files,
    new_histograms,
)
from legendsimflow.scripts import log_script_invocation
from legendsimflow.selections import (
    aoe_categories,
    aoe_flags,
    aoe_model,
    category_levels,
    fccd_by_uid,
    fccd_categories,
    has_hpge_multiplicity,
    is_hpge_hit_above_thr,
    is_spms_coincident,
    lar_categories,
    one_hit,
)


@snakemake_compatible(
    mapping={
        "evt_files": "input.evt_files",
        "simstat_part_file": "input.simstat_part_file",
        "runlist": "params.runlist",
        "pdf_file": "output.pdf_file",
        "run_dir": "output.run_dir",
        "log_file": "log[0]",
        "simid": "wildcards.simid",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the pdf tier files of the runs of a simulation, and their sum."
    )
    parser.add_argument(
        "--evt-files", nargs="+", required=True, help="evt tier files of the simulation"
    )
    parser.add_argument(
        "--simstat-part-file",
        required=True,
        help="simulation statistics partition file of the simulation",
    )
    parser.add_argument(
        "--runlist", nargs="+", required=True, help="runs of the simulation"
    )
    parser.add_argument(
        "--pdf-file", required=True, help="output pdf tier file, summed over the runs"
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="output folder of the pdf tier files of the runs",
    )
    parser.add_argument("--simid", required=True, help="simulation ID wildcard")
    parser.add_argument("--log-file", default=None, help="log file")
    parser.add_argument(
        "--simflow-config",
        "--config",
        dest="simflow_config",
        required=True,
        help="simflow config YAML path",
    )
    args = parser.parse_args()

    config = utils.init_simflow_context(args.simflow_config, workflow=None).config

    tier_pdf_settings = get_tier_settings(config, "pdf")
    tier_evt_settings = get_tier_settings(config, "evt")
    tier_hit_settings = get_tier_settings(config, "hit")

    log = ldfs.utils.build_log(
        config.metadata.simprod.config.logging.to_dict(), args.log_file
    )
    log_script_invocation(log, "tier-pdf", parser, args)

    partitions = load_dict(nersc.dvs_ro(config, args.simstat_part_file))
    runlist = sorted(set(args.runlist))
    if extra := {r for part in partitions.values() for r in part} - set(runlist):
        msg = f"runs {sorted(extra)} of the partition file are not in the run list"
        raise SimflowConfigError(msg)

    # the evt tier stores only the period and run of each event
    runid_of = {}
    for runid in runlist:
        key = parse_runid(runid)[1:3]
        if key in runid_of:
            msg = (
                f"{runid_of[key]} and {runid} differ only by the data type, "
                "their events cannot be told apart in the evt tier"
            )
            raise SimflowConfigError(msg)
        runid_of[key] = runid

    # the last job holding events of each run: its file is written after it
    last_job = {}
    for job in sorted(partitions):
        for runid in partitions[job]:
            last_job[runid] = job

    evt_files = {}
    for f in nersc.dvs_ro(config, list(args.evt_files)):
        m = re.search(r"job_\d+", Path(f).name)
        if m is None or m.group() in evt_files:
            msg = f"{f} has no job number or the same as another evt file"
            raise RuntimeError(msg)
        evt_files[m.group()] = str(f)
    if set(evt_files) != set(partitions):
        msg = (
            f"the evt files are of jobs {sorted(evt_files)}, "
            f"the partition file of jobs {sorted(partitions)}"
        )
        raise RuntimeError(msg)

    # the number of simulated primaries of each job is written by remage at
    # the root of its stp file and forwarded by the evt tier
    n_primaries = sum(
        int(lh5.read("number_of_simulated_events", f).value) for f in evt_files.values()
    )
    livetimes = {r: runinfo(config.metadata, r).livetime_in_s for r in runlist}
    primaries = primaries_per_run(livetimes, n_primaries)
    log.info(
        "... %s simulated events shared among the runs: %s", n_primaries, primaries
    )

    run_dir = Path(args.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    def run_file(runid: str) -> Path:
        name = patterns.output_tier_pdf_run_filename(
            config, simid=args.simid, runid=runid
        ).name
        return run_dir / name

    def write(runid: str, output: Struct) -> None:
        f = run_file(runid)
        lh5.write(output, "pdf", f, wo_mode="w")
        lh5.write(Scalar(primaries[runid]), "nr_sim_events", f, wo_mode="append")

    # without HPGe data (skip_hit) there is nothing to histogram
    if tier_evt_settings.skip_hit:
        for runid in runlist:
            write(runid, Struct())
    else:
        has_spms_coinc = not tier_evt_settings.skip_opt
        model = aoe_model(tier_hit_settings)
        fccd_thr = tier_pdf_settings.fccd_tag_energy_thr_kev
        if fccd_thr is not None and not tier_hit_settings.store_surface_steps:
            msg = (
                "fccd_tag_energy_thr_kev needs the steps near the n+ surface, "
                "run the hit tier with store_surface_steps"
            )
            raise ValueError(msg)

        # name to UID map written at evt-tier build time
        name_to_uid = utils.union_detector_uids(list(evt_files.values()))
        group_uids = group_detectors(tier_pdf_settings.detector_groups, name_to_uid)
        if fccd_thr is not None:
            sorted_uids, fccd_mm = fccd_by_uid(
                config.metadata,
                name_to_uid,
                tier_hit_settings.surface_depth_in_mm,
            )

        levels = category_levels(has_spms_coinc, model, fccd_thr)
        thr = tier_evt_settings.geds_energy_thr_kev

        def fill(data: ak.Array, histograms: dict, mul2_hist) -> None:
            is_hit = is_hpge_hit_above_thr(data, thr)
            energy_m2 = data.geds.energy[is_hit][
                has_hpge_multiplicity(data, 2, None, thr)
            ]
            mul2_hist.fill(ak.min(energy_m2, axis=-1), ak.max(energy_m2, axis=-1))

            # multiplicity-1 events have exactly one hit above threshold
            evt_m1 = data[has_hpge_multiplicity(data, 1, None, thr)]
            hit = one_hit(evt_m1, thr)
            energy = ak.to_numpy(hit.energy)
            uid = ak.to_numpy(hit.rawid)

            # event masks of the categories of all levels
            in_category = {}
            if has_spms_coinc:
                in_category |= lar_categories(
                    is_spms_coincident(
                        evt_m1,
                        tier_evt_settings.lar_veto_multiplicity_thr,
                        tier_evt_settings.lar_veto_energy_sum_pe_thr,
                    )
                )
            if model is not None:
                in_category |= aoe_categories(*aoe_flags(evt_m1, model, thr))
            if fccd_thr is not None:
                fccd_of_hit = fccd_mm[np.searchsorted(sorted_uids, uid)]
                in_category |= fccd_categories(evt_m1, fccd_of_hit, fccd_thr, thr)

            fill_histograms(histograms, energy, uid, in_category, group_uids)

        # read only the fields used here: the SiPM table and the surface steps
        # are large
        fields = [
            "trigger/period",
            "trigger/run",
            "geds/energy",
            "geds/rawid",
            "geds/quality/is_good_channel",
        ]
        if has_spms_coinc:
            fields += ["spms/multiplicity", "spms/energy_sum"]
        if model is not None:
            fields += ["geds/psd/is_good", f"geds/psd/{model}"]
        if fccd_thr is not None:
            fields += ["geds/surface/edep", "geds/surface/distance_to_nplus"]

        open_runs = {}
        written = set()

        def flush(done: Iterable[str]) -> None:
            for runid in done:
                log.info("... writing %s", run_file(runid))
                histograms, mul2_hist = open_runs.pop(runid)
                write(runid, format_pdf_struct(histograms, levels, mul2_hist))
                written.add(runid)

        # the evt tier writes the events of each job run by run, in the order
        # of the runs, and a run spans consecutive jobs. A run is complete when
        # a later run starts in its last job, so only the histograms of a few
        # runs are in memory at a time
        for job in sorted(evt_files):
            log.info("... reading %s", evt_files[job])
            for chunk in lh5.LH5Iterator(
                evt_files[job],
                "evt",
                field_mask=fields,
                buffer_len=tier_pdf_settings.buffer_len,
            ):
                data = chunk.view_as("ak")
                # the evt tier takes the period and run from the hit tier: they
                # are NaN for events without HPGe hits, which fill no histogram
                data = data[~np.isnan(data.trigger.run)]
                period = ak.to_numpy(data.trigger.period).astype(int)
                run = ak.to_numpy(data.trigger.run).astype(int)
                runids = {}
                for p, r in set(zip(period, run, strict=True)):
                    runid = runid_of.get((p, r))
                    if runid not in partitions[job] or runid in written:
                        msg = (
                            f"{evt_files[job]} has events of period {p} run {r}, "
                            "not in the partition file or not in the order of the runs"
                        )
                        raise RuntimeError(msg)
                    runids[runid] = (period == p) & (run == r)

                if not runids:
                    continue
                last = max(runids)
                flush([r for r in open_runs if r < min(runids) and last_job[r] == job])
                for runid in sorted(runids):
                    if runid not in open_runs:
                        open_runs[runid] = new_histograms(levels, group_uids)
                    fill(data[runids[runid]], *open_runs[runid])
                    if runid < last and last_job[runid] == job:
                        flush([runid])

            flush([r for r in open_runs if last_job[r] == job])

        # runs without events get empty histograms
        for runid in sorted(set(runlist) - written):
            histograms, mul2_hist = open_runs.pop(runid, None) or new_histograms(
                levels, group_uids
            )
            write(runid, format_pdf_struct(histograms, levels, mul2_hist))

    log.info("... summing the runs into %s", args.pdf_file)
    merge_pdf_files([run_file(r) for r in runlist], args.pdf_file)


if __name__ == "__main__":
    main()
