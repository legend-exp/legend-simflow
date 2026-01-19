# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>
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
from __future__ import annotations

from collections.abc import Iterable, Mapping


def partition_simstat(
    n_events: Mapping[str, int],
    n_events_part: Mapping[str, int],
    runlist: Iterable[str],
) -> dict[str, dict[str, list[int]]]:
    """Partition the simulation event statistics according to run livetime.

    Returns the following dictionary:

    .. code-block:: yaml

       job_000:
         l200-p03-r001-phy: [0, 300]  # interval includes its edges
         l200-p03-r002-phy: [301, 456]
       job_001:
         l200-p03-r002-phy: [0, 200]
         l200-p03-r003-phy: [201, 156]
       ...

    where the number of events of each job is partitioned in runs, such that
    the global event partitioning in `n_events_part` is respected.

    Parameters
    ----------
    n_events
        mapping of number of simulation events and simulation job.

        .. code-block:: yaml

           job_0000: 5000
           job_0001: 7000
           ...

    n_events_part
        mapping of fraction of total number of simulation events (summed over
        all jobs) per considered run, with weights equal to the run livetime
        fraction.

        .. code-block:: yaml

           l200-p03-r001-phy: 300
           l200-p03-r002-phy: 456
           ...
           l200-<...>: tot_n_events

    runlist
        list of runs in the form ``<experiment>-<period>-<run>-<datatype>``.
    """
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
            msg = (
                f"could not allocate all events in {job_id}: {job_remaining} left over"
            )
            raise RuntimeError(msg)

    # sanity: all run quotas should be exhausted
    unassigned = {r: n for r, n in remaining_run_events.items() if n != 0}
    if unassigned:
        msg = f"unassigned events left in runs: {unassigned}"
        raise RuntimeError(msg)

    return job_partitions
