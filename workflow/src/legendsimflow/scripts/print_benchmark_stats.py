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

# ruff: noqa: F821, T201

import csv
from pathlib import Path

from utils import patterns


def printline(*line):
    print("{:<52}{:>16}{:>11}{:>23}".format(*line))


printline("simid", "CPU time [ms/ev]", "evts / 1h", "jobs (1h) / 10^8 evts")
printline("-----", "----------------", "---------", "---------------------")

bdir = Path(patterns.as_ro(snakemake.config, snakemake.config["paths"]["benchmarks"]))

for simd in sorted(bdir.glob("*/*")):
    if simd.parent.name not in ("ver", "stp"):
        continue

    data = {"cpu_time": 0}
    for jobd in simd.glob("*.tsv"):
        with jobd.open(newline="") as f:
            this_data = next(iter(csv.DictReader(f, delimiter="\t")))
            data["cpu_time"] += float(this_data["cpu_time"])

    speed = (
        data["cpu_time"]
        / snakemake.config["benchmark"]["n_primaries"][simd.parent.name]
    )
    evts_1h = int(60 * 60 / speed) if speed > 0 else "..."
    njobs = int(1e8 / evts_1h) if not isinstance(evts_1h, str) else 0
    printline(
        simd.parent.name + "." + simd.name,
        "({:}s) {:.2f}".format(int(data["cpu_time"]), 1000 * speed),
        evts_1h,
        njobs,
    )
