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

import os
from pathlib import Path

from snakemake.scripts import snakemake

# TODO: generate the files from metadata inplace here
# this script is run for each data taking run
# snakemake.input[0] points to the legend-metadata clone in legend-simflow
# snakemake.output[0] is the output filename

src = (
    Path(snakemake.input[0])
    / "simprod"
    / "config"
    / "tier"
    / "evt"
    / snakemake.config["experiment"]
    / f"{snakemake.wildcards.runid}-build_evt.json"
)
dest = snakemake.output[0]

os.system(f"cp {src} {dest}")
