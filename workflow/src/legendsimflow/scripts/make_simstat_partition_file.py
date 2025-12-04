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
from legendsimflow.partitioning import partition_simstat

args = snakemake  # noqa: F821

stp_files = sorted(args.input.stp_files)
runinfo = args.params.runinfo
output_file = Path(args.output[0])
log_file = args.log[0]
metadata = args.config.metadata
runlist = sorted(args.params.runlist)

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

# get the run livetimes
run_livetimes = {run: mutils.runinfo(metadata, run).livetime_in_s for run in runlist}

# now partition the number of events according to the livetime fractions
n_events_part = {
    run: int(sum(n_events.values()) * dtime / sum(run_livetimes.values()))
    for run, dtime in run_livetimes.items()
}

# check: leftover after flooring?
leftover = sum(n_events.values()) - sum(n_events_part.values())
if leftover:
    last = next(reversed(n_events_part))
    n_events_part[last] += leftover

# create the partitioning dictionary and write it to disk
dbetto.utils.write_dict(
    partition_simstat(n_events, n_events_part, runlist), output_file
)
