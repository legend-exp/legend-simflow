# ruff: noqa: I002

# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>,
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

import argparse
import shutil

import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
from lgdo import lh5
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, utils
from legendsimflow.metadata import get_tier_settings
from legendsimflow.scripts import log_script_invocation


@snakemake_compatible(
    mapping={
        "evt_files": "input",
        "cvt_file": "output[0]",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(description="Build the cvt tier.")
    parser.add_argument(
        "--evt-files", nargs="+", required=True, help="input evt tier files"
    )
    parser.add_argument("--cvt-file", required=True, help="output cvt tier file")
    parser.add_argument("--log-file", default=None, help="log file")
    parser.add_argument(
        "--simflow-config",
        "--config",
        dest="simflow_config",
        required=True,
        help="simflow config YAML path",
    )
    args = parser.parse_args()

    config = utils.init_simflow_context(args.simflow_config, workflow=None).config

    evt_files = nersc.dvs_ro(config, list(args.evt_files))
    cvt_file, move2cfs = nersc.make_on_scratch(config, args.cvt_file)
    log_file = args.log_file

    tier_cvt_settings = get_tier_settings(config, "cvt")
    buffer_len = tier_cvt_settings.buffer_len

    log = ldfs.utils.build_log(config.metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "tier-cvt", parser, args)
    log.info("starting cvt merge")

    if len(evt_files) == 1:
        shutil.copy(evt_files[0], cvt_file)
    else:
        # detector_uids must be identical across all evt inputs; merging would
        # silently mask mismatches, so assert before copying.
        ref = lh5.read("detector_uids", evt_files[0])
        ref_map = {name: int(ref[name].value) for name in ref}
        for f in evt_files[1:]:
            other = lh5.read("detector_uids", f)
            other_map = {name: int(other[name].value) for name in other}
            if other_map != ref_map:
                msg = (
                    f"detector_uids in {f} disagrees with {evt_files[0]}; "
                    "cvt merge requires all evt inputs to share the same "
                    "name->rawid mapping"
                )
                raise ValueError(msg)
        lh5.write(ref, "detector_uids", cvt_file, wo_mode="write_safe")

        for table in lh5.ls(evt_files[0]):
            if table == "detector_uids":
                continue
            for chunk in lh5.LH5Iterator(evt_files, table, buffer_len=buffer_len):
                lh5.write(chunk, table, cvt_file, wo_mode="append")

    move2cfs()


if __name__ == "__main__":
    main()
