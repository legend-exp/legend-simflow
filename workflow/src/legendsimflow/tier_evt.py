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


def smk_get_evt_window(wildcards, input):
    # open file with run livetime partitioning
    with Path(input.run_part_file[0]).open() as f:
        runpart = json.load(f, object_pairs_hook=OrderedDict)

    runs = list(runpart.keys())
    weights = list(runpart.values())

    # compute start event and number of events for this run
    start_event = sum(weights[: runs.index(wildcards.runid)])
    n_events = runpart[wildcards.runid]

    return (int(start_event), int(n_events))
