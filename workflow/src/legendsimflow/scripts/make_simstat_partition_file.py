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

import re
from pathlib import Path

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
from lgdo import lh5

from legendsimflow import metadata as mutils

stp_files = sorted(snakemake.input.stp_files)  # noqa: F821
runinfo = snakemake.params.runinfo  # noqa: F821
output_file = Path(snakemake.output[0])  # noqa: F821
log_file = snakemake.log[0]  # noqa: F821
metadata = snakemake.config.metadata  # noqa: F821
runlist = sorted(snakemake.params.runlist)  # noqa: F821

# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)

# get full simulation event statistics
n_events = {}
for stp_file in stp_files:
    m = re.search(r"job_\d+", stp_file)
    if m is None:
        msg = f"{stp_file} does not contain the job number segment"
        raise RuntimeError(msg)

    n_events[m.group()] = lh5.utils.read_n_rows("/tcm", stp_file)

tot_n_events = sum(n_events.values())

# get the run livetimes
run_livetimes = {run: mutils.runinfo(metadata, run).livetime_in_s for run in runlist}
tot_livetime = sum(run_livetimes.values())

# now partition the number of events according to the livetime fractions
n_events_part = {
    run: int(tot_n_events * dtime / tot_livetime)
    for run, dtime in run_livetimes.items()
}

# leftover after flooring?
leftover = tot_n_events - sum(n_events_part.values())
if leftover:
    last = next(reversed(n_events_part))
    n_events_part[last] += leftover

# this is now in the format:
#
# l200-p03-r001-phy: 300
# l200-p03-r002-phy: 456
# ...
# l200-<...>: tot_n_events

# we want to build the following dictionary:
#
# job_000:
#   l200-p03-r001-phy: [0, 300] # interval includes its edges
#   l200-p03-r002-phy: [301, 456]
# job_001:
#   l200-p03-r002-phy: [0, 200]
#   l200-p03-r003-phy: [201, 156]
# ...
#
# where the number of events in each file is partitioned in runs, such that the
# global partitioning in n_events_part is respected.

# remaining events per run to allocate across jobs
remaining_run_events = dict(n_events_part)

# job -> run -> [start, end]
job_partitions: dict[str, dict[str, list[int]]] = {}

for job_id, job_nevents in n_events.items():
    job_remaining = job_nevents
    start = 0
    job_partitions[job_id] = {}

    for run in runlist:
        # if we have consumed all events from this job, let's move to the next
        if job_remaining <= 0:
            break

        # how many events do we still need to consume from this partition?
        run_remaining = remaining_run_events[run]
        if run_remaining <= 0:
            continue

        # consume all events needed for this run or all remaining events from the job
        take = min(job_remaining, run_remaining)
        # register the chunk
        job_partitions[job_id][run] = [start, start + take - 1]  # [start, end]
        # prepare for the next iteration
        start += take
        job_remaining -= take
        remaining_run_events[run] -= take

    if job_remaining != 0:
        msg = f"could not allocate all events in {job_id}: {job_remaining} left over"
        raise RuntimeError(msg)

# sanity: all run quotas should be exhausted
unassigned = {r: n for r, n in remaining_run_events.items() if n != 0}
if unassigned:
    msg = f"unassigned events left in runs: {unassigned}"
    raise RuntimeError(msg)

dbetto.utils.write_dict(job_partitions, output_file)
