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

import hashlib
import inspect
import json
import logging
import os
from datetime import datetime
from pathlib import Path

import dbetto
import legenddataflowscripts as ldfs
from dbetto import AttrsDict
from legendmeta import LegendMetadata
from numpy.typing import ArrayLike
from reboost.hpge.psd import _current_pulse_model as current_pulse_model

from . import SimflowConfig, nersc

log = logging.getLogger(__name__)


def _merge_defaults(user: dict, default: dict) -> dict:
    # merge default into user without overwriting user values
    result = dict(default)
    for k, v in user.items():
        if k in result and isinstance(result[k], dict) and isinstance(v, dict):
            result[k] = _merge_defaults(v, result[k])
        else:
            result[k] = v
    return result


def init_simflow_context(raw_config: dict, workflow) -> AttrsDict:
    """Pre-process and sanitize the Simflow configuration.

    - set default configuration fields;
    - substitute ``$_`` and environment variables;
    - convert to :class:`~dbetto.attrsdict.AttrsDict`;
    - cast filesystem paths to :class:`pathlib.Path`;
    - clone and configure `legend-metadata`;
    - attach a :class:`legendmeta.LegendMetadata` instance to the Simflow
      configuration;
    - export important environment variables.

    Returns a dictionary with useful objects to be used in the Simflow
    Snakefiles (i.e. the "context").
    """
    if not raw_config:
        msg = "you must set a config file with --configfile"
        raise RuntimeError(msg)

    raw_config = _merge_defaults(
        raw_config,
        {"benchmark": {"enabled": False}, "nersc": {"dvs_ro": False, "scratch": False}},
    )

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

    if "l200data" in config.paths:
        config["paths"]["l200data"] = nersc.dvs_ro(config, config.paths.l200data)

    # NOTE: this will attempt a clone of legend-metadata, if the directory does not exist
    # NOTE: don't use lazy=True, we need a fully functional TextDB
    # NOTE: read only path on NERSC, we are not going to modify the db
    metadata = LegendMetadata(nersc.dvs_ro(config, config.paths.metadata))

    if "legend_metadata_version" in config:
        metadata.checkout(config.legend_metadata_version)

    config["metadata"] = metadata

    # make sure all simflow plots are made with a consistent style
    # I have verified only that this variable is visible in scripts (not shell directives)
    os.environ["MPLCONFIGDIR"] = f"{workflow.basedir}/src/legendsimflow"

    proctime = (
        "benchmark"
        if config.benchmark.enabled
        else datetime.now().strftime("%Y%m%dT%H%M%SZ")
    )

    return AttrsDict(
        {
            "config": config,
            "basedir": workflow.basedir,
            "proctime": proctime,
        }
    )


def setup_logdir_link(config: SimflowConfig, proctime):
    """Set up the timestamp-tagged directory for the workflow log files."""
    logdir = Path(config.paths.log)
    logdir.mkdir(parents=True, exist_ok=True)

    # create a handy link to access latest log directory
    link = logdir / "latest"
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(proctime, target_is_directory=True)


def lookup_dataflow_config(l200data: Path | str) -> AttrsDict:
    """Find the paths to the data inputs.

    Parameters
    ----------
    l200data
        The path to the L200 data production cycle.
    Returns
    -------
    the dataflow configuration file as a dictionary.
    """
    if not isinstance(l200data, Path):
        l200data = Path(l200data)

    # look for the dataflow config file
    df_cfgs = [p for p in l200data.glob("*config.*") if not p.name.startswith(".")]

    if len(df_cfgs) > 1:
        msg = f"found multiple configuration files in {l200data}, this cannot be!"
        raise RuntimeError(msg)

    if len(df_cfgs) == 0:
        msg = f"did not find any valid config file in {l200data}, this cannot be!"
        raise RuntimeError(msg)

    msg = f"found dataflow configuration file: {df_cfgs[0]}"
    log.debug(msg)

    return dbetto.utils.load_dict(df_cfgs[0])


def get_hit_tier_name(l200data: str) -> str:
    """Extract the name of the hit tier for this production cycle.

    If the `pht` tier is present this is used else the `hit` tier is used.

    Parameters
    ----------
    l200data
        Path to the production cycle of l200 data.
    """

    dataflow_config = lookup_dataflow_config(l200data)

    df_cfg = (
        dataflow_config["setups"]["l200"]["paths"]
        if ("setups" in dataflow_config)
        else dataflow_config["paths"]
    )

    # first check if pht exists
    has_pht = ("tier_pht" in df_cfg) and Path(
        df_cfg["tier_pht"].replace("$_", str(l200data))
    ).exists()
    has_hit = ("tier_hit" in df_cfg) and Path(
        df_cfg["tier_hit"].replace("$_", str(l200data))
    ).exists()

    if has_pht:
        return "pht"
    if has_hit:
        return "hit"
    msg = f"The l200data {l200data} does not contain a valid pht or hit tier"
    raise RuntimeError(msg)


def hash_dict(d: dict | AttrsDict) -> str:
    """Compute the hash of a Python dict."""
    if isinstance(d, AttrsDict):
        d = d.to_dict()

    return json.dumps(d, sort_keys=True)

    # NOTE: alternatively, return sha256 (shorter string but bad for diffs)
    # return hashlib.sha256(s.encode()).hexdigest()


def string_to_remage_seed(s: str, seed: int = 0) -> int:
    """Convert a string to a positive 32-bit integer.

    This is good to use as remage seed.

    Parameters
    ----------
    s
        the string to encode.
    seed
        optional seed.
    """
    seed_bytes = seed.to_bytes(8, byteorder="big", signed=False)  # 0..2^64-1
    h = hashlib.sha256(seed_bytes + s.encode("utf-8")).digest()
    return int.from_bytes(h[:4], "big", signed=False) & 0x7FFFFFFF


def _curve_fit_popt_to_dict(popt: ArrayLike) -> dict:
    """Get the :func:`scipy.curve_fit` parameter results as a dictionary"""
    params = list(inspect.signature(current_pulse_model).parameters)
    param_names = params[1:]

    popt_dict = dict(zip(param_names, popt, strict=True))

    for key, value in popt_dict.items():
        popt_dict[key] = float(f"{value:.3g}")

    return popt_dict
