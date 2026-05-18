# ruff: noqa: I002

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

import argparse

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import metadata as mutils
from legendsimflow import utils
from legendsimflow.scripts import log_script_invocation


@snakemake_compatible(
    mapping={
        "runid": "wildcards.runid",
        "hpge_detector": "wildcards.hpge_detector",
        "pars_file": "output.pars_file",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract the HPGe electronics model for a LEGEND run."
    )
    parser.add_argument(
        "--runid",
        required=True,
        help="LEGEND run identifier (e.g. l200-p03-r000-phy)",
    )
    parser.add_argument(
        "--hpge-detector",
        required=True,
        dest="hpge_detector",
        help="HPGe detector name",
    )
    parser.add_argument(
        "--pars-file",
        required=True,
        help="output YAML file for the electronics model parameters",
    )

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
    metadata = config.metadata

    log = ldfs.utils.build_log(metadata.simprod.config.logging, args.log_file)
    log_script_invocation(log, "extract-hpge-elecmod", parser, args)

    runid = args.runid
    hpge = args.hpge_detector
    pars_file = args.pars_file

    # check for metadata-driven defaults first; if present, bypass the
    # waveform fitting entirely (enables LEGEND-1000 simulations without
    # l200data)
    raw_elecmod = mutils.simpars(
        metadata, "geds.elecmod", runid, config.experiment, default=None
    )
    elecmod_default = (
        raw_elecmod.get("default", None) if raw_elecmod is not None else None
    )

    if elecmod_default is not None:
        log.info("... using elecmod metadata defaults for %s in %s", hpge, runid)
        entry = raw_elecmod.get(hpge, elecmod_default)
        dbetto.utils.write_dict(entry.to_dict(), pars_file)
        return

    msg = "Currently only metadata driven electronics model is supported."
    raise NotImplementedError(msg)


if __name__ == "__main__":
    main()
