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

import math
import re
from datetime import timedelta
from statistics import mean

from ruamel.yaml import YAML

from legendsimflow import nersc


def round_down_2sf_5(x):
    k = math.floor(math.log10(x))
    m = x / 10**k

    m_rd = round(m / 0.5) * 0.5
    return m_rd * 10**k


def printline(*line):
    print("{:<52}{:>16}{:>27}{:>11}{:>12}{:>23}{:>12}".format(*line))


args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

speed_pattern = re.compile(
    r"^.*Stats: average event processing time was\s+"
    r"([0-9]+(?:\.[0-9]+)?)\s+seconds/event\s+=\s+"
    r"([0-9]+(?:\.[0-9]+)?)\s+events/second\s*$",
    re.MULTILINE,
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
    "...rounded",
    "jobs (1h) / 10^8 evts",
    "...rounded",
)
printline(
    "-----",
    "-------------",
    "-------------------------",
    "---------",
    "----------",
    "---------------------",
    "----------",
)

updates: dict[str, tuple[int, int]] = {}

for simd in sorted(logdir.glob("*/*")):
    # this code works only for remage output
    if simd.parent.name != "stp":
        continue

    speed = "..."
    runtime = "..."
    evts_1h = "..."
    njobs = "..."
    evts_1h_round = "..."
    njobs_round = "..."

    # pool the per-thread samples of all benchmark jobs of this simid
    thread_speeds: list[float] = []
    runtimes: list[float] = []

    for jobd in simd.glob("*.log"):
        with jobd.open("r", encoding="utf-8") as f:
            data = f.read()

        # extract events/sec for each thread
        job_speeds = [float(m.group(2)) for m in speed_pattern.finditer(data)]

        # get the runtime of each thread
        job_runtimes = [
            timedelta(
                days=int(d), hours=int(h), minutes=int(mi), seconds=int(s)
            ).total_seconds()
            for d, h, mi, s in time_pattern.findall(data)
        ]

        # simulation might have crashed or still be running
        if not job_speeds or not job_runtimes:
            continue

        thread_speeds += job_speeds
        runtimes += job_runtimes

    if thread_speeds and runtimes:
        runtime = mean(runtimes)
        speed = mean(thread_speeds)

        evts_1h = int(speed * 3600)
        evts_1h_round = int(round_down_2sf_5(speed * 3600))
        njobs = int(1e8 / evts_1h)
        njobs_round = int(1e8 / evts_1h_round)

    printline(
        simd.parent.name + "." + simd.name,
        ("!!! " if isinstance(runtime, float) and runtime < 10 else "")
        + (f"{runtime:.1f}" if isinstance(runtime, float) else runtime),
        f"{speed:.2f}" if isinstance(speed, float) else speed,
        evts_1h,
        evts_1h_round,
        njobs,
        njobs_round,
    )

    if isinstance(evts_1h_round, int):
        updates[simd.name] = evts_1h_round

# generate updated simconfig.yaml
simconfig_path = (
    args.config.paths.config / "tier/stp" / args.config.experiment / "simconfig.yaml"
)

ryaml = YAML()
ryaml.preserve_quotes = True
with simconfig_path.open(encoding="utf-8") as f:
    simconfig = ryaml.load(f)

target_evts = int(1e7)
for simid, evts_per_job in updates.items():
    simconfig[simid]["primaries_per_job"] = evts_per_job
    simconfig[simid]["number_of_jobs"] = max(1, math.ceil(target_evts / evts_per_job))

out = args.config.paths.benchmarks / "generated-simconfig.yaml"
out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", encoding="utf-8") as f:
    ryaml.dump(simconfig, f)

print(f"\nGenerated simconfig (1h jobs, >={target_evts:.0e} total events): {out}")
