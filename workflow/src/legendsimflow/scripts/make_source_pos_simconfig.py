# Copyright (C) 2026 Toby Dixon <toby.dixon.23@ucl.ac.uk>,
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

"""
Update the simconfig for the source position analysis.

Based on a template configuration (in input_path) which is expanded
out in the output_path based on the grid defined by z and phi.

Usage:
    python make_source_pos_simconfig.py [options]
"""

from __future__ import annotations

import argparse
import copy
import logging
import shutil
from pathlib import Path

import dbetto

log = logging.getLogger(__name__)


PRIMARIES = 4_000_000
JOBS = 5
ONLY_TH = True


def replace_position(config: dict, dz: float, dphi: float):
    # validate required geometry/SIS structure
    geom_extra = config.get("geom_config_extra")
    if not isinstance(geom_extra, dict):
        msg = (
            "Template config must contain a 'geom_config_extra' mapping with a 'sis' "
            "subtree defining exactly one SIS for the source position analysis."
        )
        raise ValueError(msg)

    if "sis" not in geom_extra or not isinstance(geom_extra["sis"], dict):
        msg = (
            "Template config must contain 'geom_config_extra[\"sis\"]' as a mapping "
            "defining exactly one SIS for the source position analysis."
        )
        raise ValueError(msg)

    sis_config = geom_extra["sis"]
    if not sis_config:
        msg = (
            "Template config must define at least one SIS under "
            "'geom_config_extra[\"sis\"]' for the source position analysis."
        )
        raise ValueError(msg)

    # get the sis
    if len(list(sis_config.keys())) != 1:
        msg = "Can only have source in one SIS!"
        raise ValueError(msg)

    sis = next(iter(sis_config))

    config["primaries_per_job"] = PRIMARIES
    config["number_of_jobs"] = JOBS

    sis_config[sis]["offset"] = dz
    sis_config[sis]["phi_offset"] = dphi

    return config


def main():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(description="Process arguments.")

    # positional string argument
    parser.add_argument("--input-path", "-i", type=str, help="Input file path")
    parser.add_argument("--output-path", "-o", type=str, help="Output file path")

    # list arguments
    parser.add_argument(
        "--phi", type=float, nargs="+", required=True, help="List of phi values"
    )
    parser.add_argument(
        "--z", type=float, nargs="+", required=True, help="List of z values"
    )

    args = parser.parse_args()

    z = args.z
    phi = args.phi

    configuration = Path(args.input_path)

    log.info("... copying input %s to %s", args.input_path, args.output_path)
    # create the output

    if args.input_path == args.output_path:
        msg = "Input and output path cannot be the same!"
        raise ValueError(msg)

    # also replace the hit config
    hit_input = args.input_path.replace("stp", "hit")
    hit_output = args.output_path.replace("stp", "hit")

    for indir, outdir in zip(
        [hit_input, args.input_path], [hit_output, args.output_path], strict=True
    ):
        if Path(outdir).exists():
            shutil.rmtree(Path(outdir))

        shutil.copytree(indir, outdir)

    log.info("... editing simconfig based on z: %s and phi: %s", z, phi)

    # now edit the simconfig
    input_simconfig = dbetto.utils.load_dict(configuration / "simconfig.yaml")
    output_simconfig = {}

    input_hit_simconfig = dbetto.utils.load_dict(Path(hit_input) / "simconfig.yaml")
    output_hit_simconfig = {}

    for name, info in input_simconfig.items():
        # only th for now
        if "Pb214_to_Po214" in name and ONLY_TH:
            continue

        for ztemp in z:
            for phi_temp in phi:
                new_block = replace_position(
                    copy.deepcopy(info), dz=ztemp, dphi=phi_temp
                )
                output_simconfig[
                    f"{name}_z_{int(ztemp)}_mm_phi_{int(phi_temp)}_deg"
                ] = new_block
                output_hit_simconfig[
                    f"{name}_z_{int(ztemp)}_mm_phi_{int(phi_temp)}_deg"
                ] = copy.deepcopy(input_hit_simconfig[name])

    log.info("... saving the result in %s/simconfig.yaml", args.output_path)

    dbetto.utils.write_dict(output_simconfig, Path(args.output_path) / "simconfig.yaml")
    dbetto.utils.write_dict(output_hit_simconfig, Path(hit_output) / "simconfig.yaml")


if __name__ == "__main__":
    main()
