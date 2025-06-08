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

import json
from pathlib import Path

import uproot
from utils import utils

with Path(snakemake.input.runinfo).open() as f:
    rinfo = json.load(f)

# retrieve run livetimes
runs = utils.get_some_list(list(set(snakemake.config["runlist"])))
spec = [r.split("-") for r in runs]
livetimes = [rinfo[p[1]][p[2]][p[3]]["livetime_in_s"] for p in spec]

# get total number of mc events from hit files
print(
    "INFO: computing total number of MC events for simulation "
    + snakemake.wildcards.simid
)
file_evts = uproot.num_entries(
    [f"{file}:simTree" for file in snakemake.params.ro_hit_files]
)
tot_events = 0
for file in file_evts:
    tot_events += file[-1]

# write json file with weights for each run
with Path(snakemake.output[0]).open("w") as f:
    json.dump(
        dict(zip(runs, [int(tot_events * t / sum(livetimes)) for t in livetimes])),
        f,
        indent=2,
    )
