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

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import legenddataflowscripts as ldfs
from dbetto import AttrsDict
from legendmeta import LegendMetadata

from . import SimflowConfig


def init_simflow_context(raw_config: dict, workflow) -> AttrsDict:
    if not raw_config:
        msg = "you must set a config file with --configfile"
        raise RuntimeError(msg)

    raw_config.setdefault("benchmark", {"enabled": False})
    ldfs.workflow.utils.subst_vars_in_snakemake_config(workflow, raw_config)
    config = AttrsDict(raw_config)

    # convert all strings in the "paths" block to pathlib.Path
    def _make_path(d):
        for k, v in d.items():
            if isinstance(v, str):
                d[k] = Path(v)
            else:
                d[k] = _make_path(v)
        return d

    config["paths"] = _make_path(config.paths)

    # NOTE: this will attempt a clone of legend-metadata, if the directory does not exist
    # NOTE: don't use lazy=True, we need a fully functional TextDB
    metadata = LegendMetadata(config.paths.metadata)
    if "legend_metadata_version" in config:
        metadata.checkout(config.legend_metadata_version)
    config["metadata"] = metadata

    return AttrsDict(
        {
            "config": config,
            "basedir": workflow.basedir,
            "proctime": datetime.now().strftime("%Y%m%dT%H%M%SZ"),
        }
    )


def setup_logdir_link(config: SimflowConfig, proctime):
    logdir = Path(config.paths.log)
    logdir.mkdir(parents=True, exist_ok=True)

    # create a handy link to access latest log directory
    link = logdir / "latest"
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(proctime, target_is_directory=True)
