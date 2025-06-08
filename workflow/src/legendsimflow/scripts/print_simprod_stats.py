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
import json
from datetime import timedelta
from pathlib import Path

from snakemake.scripts import snakemake
from utils import patterns


def printline(*line):
    print("{:<50} {:>20} {:>6} {:>14} {:>10}".format(*line))


printline("     ", "wall time [s]", "    ", "wall time [s]", "         ")
printline("simid", " (cumulative)", "jobs", "    (per job)", "primaries")
printline("-----", "-------------", "----", "-------------", "---------")

bdir = Path(patterns.as_ro(snakemake.config, snakemake.config["paths"]["benchmarks"]))

tot_wall_time = 0
for simd in sorted(bdir.glob("*/*")):
    njobs = 0
    data = {"wall_time": 0}
    for jobd in simd.glob("*.tsv"):
        njobs += 1
        with jobd.open(newline="") as f:
            this_data = next(iter(csv.DictReader(f, delimiter="\t")))
            data["wall_time"] += float(this_data["s"])
    tot_wall_time += data["wall_time"]

    if njobs == 0:
        continue

    tier = simd.parent.name if simd.parent.name in ("ver", "stp") else "stp"

    tdir = patterns.template_macro_dir(snakemake.config, tier=tier)
    with (tdir / "simconfig.json").open() as f:
        config = json.load(f)[simd.name]

    nprim = config["number_of_primaries"]

    printline(
        simd.parent.name + "." + simd.name,
        str(timedelta(seconds=int(data["wall_time"]))),
        njobs,
        str(timedelta(seconds=int(data["wall_time"] / njobs))),
        f"{nprim:.2E}",
    )

print("\nTotal wall time:", str(timedelta(seconds=int(tot_wall_time))))
