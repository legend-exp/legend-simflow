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
from pathlib import Path

import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
from lgdo import Scalar, Struct
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, utils
from legendsimflow.metadata import get_tier_settings
from legendsimflow.scripts import log_script_invocation


def union_detector_uids(evt_files: list[str | Path]) -> dict[str, int]:
    """Build the union ``name -> uid`` mapping over a set of evt files.

    Different jobs in the same simid can produce evt files with different
    detector subsets when low-rate decays leave some detectors with zero hits.
    Their ``detector_uids`` structs then disagree on which detector names are
    present, even though the underlying ``name -> uid`` mapping is consistent.
    This helper unions the entries; any genuine collision (a name mapped to
    different uids, or a uid mapped to different names, across files) is
    treated as a real inconsistency and raises.
    """
    union: dict[str, int] = {}
    uid_to_name: dict[int, str] = {}
    for f in evt_files:
        per_file = lh5.read("detector_uids", f)
        for name in per_file:
            uid = int(per_file[name].value)
            if name in union and union[name] != uid:
                msg = (
                    f"detector_uids inconsistency: '{name}' maps to {union[name]} "
                    f"in earlier evt files but to {uid} in {f}"
                )
                raise ValueError(msg)
            if uid in uid_to_name and uid_to_name[uid] != name:
                msg = (
                    f"detector_uids inconsistency: uid {uid} maps to "
                    f"'{uid_to_name[uid]}' in earlier evt files but to '{name}' in {f}"
                )
                raise ValueError(msg)
            union[name] = uid
            uid_to_name[uid] = name
    return union


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
        # the single evt file already carries the root-level
        # number_of_simulated_events scalar, which is copied along verbatim.
        shutil.copy(evt_files[0], cvt_file)
    else:
        # detector_uids may differ across jobs because low-rate decays can leave
        # some detectors with zero hits in some jobs (and thus absent from
        # detector_uids). We union the mappings; a real (name, uid) collision is
        # raised by union_detector_uids.
        merged_uids = union_detector_uids(list(evt_files))
        merged = Struct({name: Scalar(uid) for name, uid in merged_uids.items()})
        lh5.write(merged, "detector_uids", cvt_file, wo_mode="write_safe")

        for table in lh5.ls(evt_files[0]):
            # number_of_simulated_events is a root-level scalar summed below,
            # not a table
            if table in ("detector_uids", "number_of_simulated_events"):
                continue
            for chunk in lh5.LH5Iterator(evt_files, table, buffer_len=buffer_len):
                lh5.write(chunk, table, cvt_file, wo_mode="append")

        # total number of simulated primary events for this simid is the sum of
        # the per-job counts forwarded from the stp files by the evt tier.
        total_n_events = sum(
            int(lh5.read("number_of_simulated_events", f).value) for f in evt_files
        )
        lh5.write(
            Scalar(total_n_events),
            "number_of_simulated_events",
            cvt_file,
            wo_mode="append",
        )

    move2cfs()


if __name__ == "__main__":
    main()
