# ruff: noqa: I002, T201

# Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
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
from datetime import timedelta
from statistics import mean

from legendsimflow import nersc


def printline(*line):
    print("{:<52}{:>16}{:>27}{:>11}{:>23}".format(*line))


args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

speed_pattern = re.compile(
    r"^.*Stats: average event processing time was\s+"
    r"([0-9]+(?:\.[0-9]+)?)\s+seconds/event\s+=\s+"
    r"([0-9]+(?:\.[0-9]+)?)\s+events/second\s*$",
    re.MULTILINE,
)

nev_pattern = re.compile(
    r"^.*Run nr\. \d+ completed\. (\d+) events simulated\.", re.MULTILINE
)

time_pattern = re.compile(
    r"^.*Stats: run time was (\d+) days, (\d+) hours, (\d+) minutes and (\d+) seconds$",
    re.MULTILINE,
)

# have a look at the latest run
logdir = (nersc.dvs_ro(args.config, args.config.paths.log) / "benchmark").resolve()

if not logdir.is_dir():
    msg = "no benchmark run available!"
    raise RuntimeError(msg)

printline(
    "simid",
    "runtime [sec]",
    "speed (hot loop) [ev/sec]",
    "evts / 1h",
    "jobs (1h) / 10^8 evts",
)
printline(
    "-----",
    "-------------",
    "-------------------------",
    "---------",
    "---------------------",
)

for simd in sorted(logdir.glob("*/*")):
    # this code works only for remage output
    if simd.parent.name != "stp":
        continue

    speed = 0
    runtime = 0
    for jobd in simd.glob("*.log"):
        with jobd.open("r", encoding="utf-8") as f:
            # read the full file in memory (assuming it can't be huge)
            data = f.read()

            # extract events/sec for each thread
            time = [
                float(m.group(2)) for m in speed_pattern.finditer(data) if m is not None
            ]

            # simulations might have crashed or still running
            if time == []:
                runtime = "..."
                speed = "..."

            # get the number of simulated events for each thread (it's always the same)
            nev = int(nev_pattern.search(data).group(1))

            # get the runtime of each thread
            runtimes = [
                timedelta(
                    days=int(d), hours=int(h), minutes=int(mi), seconds=int(s)
                ).total_seconds()
                for d, h, mi, s in time_pattern.findall(data)
            ]

            runtime = mean(runtimes)
            speed += mean(time)

    evts_1h = int(speed * 60 * 60) if speed > 0 else "..."
    njobs = int(1e8 / evts_1h) if not isinstance(evts_1h, str) else 0
    printline(
        simd.parent.name + "." + simd.name,
        ("!!! " if runtime < 10 else "") + f"{runtime:.1f}",
        f"{speed:.2f}",
        evts_1h,
        njobs,
    )
