# ruff: noqa: I002, T201

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

import re
from statistics import mean

import dbetto

from legendsimflow import nersc, optimise

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

MAX_RUNTIME_MINS = args.params.max_runtime_mins
TARGET_TOTAL_PRIMARIES = args.params.target_total_primaries
N_SIG_FIGS = args.params.get("n_sig_figs", 2)

speed_pattern = re.compile(
    r"^.*Stats: average event processing time was\s+"
    r"([0-9]+(?:\.[0-9]+)?)\s+seconds/event\s+=\s+"
    r"([0-9]+(?:\.[0-9]+)?)\s+events/second\s*$",
    re.MULTILINE,
)

logdir = (nersc.dvs_ro(args.config, args.config.paths.log) / "benchmark").resolve()

if not logdir.is_dir():
    msg = "no benchmark run available!"
    raise RuntimeError(msg)


def printline(*line):
    print("{:<52}{:>25}{:>12}{:>25}{:>20}".format(*line))


printline(
    "simid",
    "primaries_per_job",
    "number_of_jobs",
    "total_primaries",
    "est. job time [min]",
)
printline(
    "-----",
    "-----------------",
    "--------------",
    "---------------",
    "-------------------",
)

settings = {}

for simd in sorted(logdir.glob("*/*")):
    # only process remage (stp tier) output
    if simd.parent.name != "stp":
        continue

    speed = 0.0
    valid = True
    for jobd in simd.glob("*.log"):
        with jobd.open("r", encoding="utf-8") as f:
            data = f.read()

        thread_speeds = [float(m.group(2)) for m in speed_pattern.finditer(data)]

        if not thread_speeds:
            valid = False
            break

        speed += mean(thread_speeds)

    if not valid or speed <= 0:
        printline(
            "stp." + simd.name,
            "N/A",
            "N/A",
            "N/A",
            "N/A",
        )
        continue

    ppj, njobs = optimise.compute_sim_settings(
        speed=speed,
        max_runtime_mins=MAX_RUNTIME_MINS,
        target_total_primaries=TARGET_TOTAL_PRIMARIES,
        n_sig_figs=N_SIG_FIGS,
    )

    total_primaries = ppj * njobs
    est_runtime_mins = ppj / speed / 60

    printline(
        "stp." + simd.name,
        f"{ppj:,}",
        str(njobs),
        f"{total_primaries:,}",
        f"{est_runtime_mins:.1f}",
    )

    settings[simd.name] = {
        "primaries_per_job": ppj,
        "number_of_jobs": njobs,
    }

print(
    f"\nConstraints: max_runtime={MAX_RUNTIME_MINS} min,"
    f" target_total_primaries={TARGET_TOTAL_PRIMARIES:.2E},"
    f" n_sig_figs={N_SIG_FIGS}"
)
print(f"Settings written to: {args.output[0]}")

dbetto.utils.write_dict(settings, args.output[0])
