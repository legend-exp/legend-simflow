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
from collections import OrderedDict
from pathlib import Path

import uproot

from . import patterns, aggregate


def smk_get_evt_window(config, wildcards, input):
    # open file with run livetime partitioning
    with Path(input.run_part_file[0]).open() as f:
        runpart = json.load(f, object_pairs_hook=OrderedDict)

    # get total number of mc events
    file_evts = uproot.num_entries([f"{file}:simTree" for file in input.hit_files])

    tot_events = 0
    for file in file_evts:
        tot_events += file[-1]

    runs = list(runpart.keys())
    weights = [tot_events * v for v in runpart.values()]

    # compute start event and number of events for this run
    start_event = sum(weights[: runs.index(wildcards.runid)])
    n_events = tot_events * runpart[wildcards.runid]

    return (int(start_event), int(n_events))