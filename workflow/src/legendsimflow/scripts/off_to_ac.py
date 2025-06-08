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

# NOTE: this script is not ready to be integrated in Snakemake, since it only
# works on the LNGS cluster! At the end of this file there's an example of how
# it should look like

# ruff: noqa: F821, T201


import json
import os
import shutil
from pathlib import Path

path = "/global/homes/t/tdixon/LEGEND/legend-simflow-config/"
outpath = "/global/homes/t/tdixon/LEGEND/legend-simflow-config-no-ac/"

sub_path = "tier/evt/l200a/"
try:
    shutil.rmtree(outpath)
    shutil.copytree(path, outpath)
except shutil.Error as e:
    print(f"Error: {e}")
except OSError as e:
    print(f"Error: {e}")

files = os.listdir(outpath + sub_path)  # noqa: PTH208
for f in files:
    if "l200" not in f:
        continue

    with Path.open(outpath + sub_path + f) as file:
        cfg = json.load(file)
    for c in cfg:
        if cfg[c]["usability"] == "off":
            cfg[c]["usability"] = "ac"
            cfg[c]["energy"] = {"sig0": 1, "sig1": 0, "sig2": 0}

    with Path.open(outpath + sub_path + f, "w") as file:
        json.dump(cfg, file, indent=2)
