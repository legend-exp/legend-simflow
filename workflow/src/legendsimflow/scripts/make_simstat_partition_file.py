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
from pathlib import Path

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import metadata as mutils
from legendsimflow import nersc, utils
from legendsimflow.partitioning import partition_simstat
from legendsimflow.scripts import log_script_invocation


@snakemake_compatible(
    mapping={
        "stp_files": "input.stp_files",
        "runlist": "params.runlist",
        "output_file": "output[0]",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the simulation statistics partition file."
    )
    parser.add_argument(
        "--stp-files", nargs="+", required=True, help="input stp tier files"
    )
    parser.add_argument(
        "--runlist", nargs="+", required=True, help="list of run IDs to partition over"
    )
    parser.add_argument(
        "--output-file", required=True, help="output partition YAML file"
    )

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

    stp_files = nersc.dvs_ro(config, sorted(args.stp_files))
    output_file = Path(args.output_file)
    log_file = args.log_file
    metadata = config.metadata
    runlist = sorted(args.runlist)

    # setup logging
    log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "make-simstat-partition-file", parser, args)

    # get full simulation event statistics
    n_events = {}
    for stp_file in stp_files:
        m = re.search(r"job_\d+", stp_file)
        if m is None:
            msg = f"{stp_file} does not contain the job number segment"
            raise RuntimeError(msg)

        n_events[m.group()] = lh5.read_n_rows("/tcm", stp_file)

    log.debug("n_events per job: %s", n_events)

    # get the run livetimes
    run_livetimes = {
        run: mutils.runinfo(metadata, run).livetime_in_s for run in runlist
    }

    # now partition the number of events according to the livetime fractions
    total_events = sum(n_events.values())
    total_livetime = sum(run_livetimes.values())
    n_events_part = {
        run: int(total_events * dtime / total_livetime)
        for run, dtime in run_livetimes.items()
    }

    # check: leftover after flooring?
    leftover = total_events - sum(n_events_part.values())
    if leftover:
        last = next(reversed(n_events_part))
        n_events_part[last] += leftover

    # create the partitioning dictionary and write it to disk
    dbetto.utils.write_dict(
        partition_simstat(n_events, n_events_part, runlist), output_file
    )


if __name__ == "__main__":
    main()
