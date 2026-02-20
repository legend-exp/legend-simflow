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
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path

import awkward as ak
import dbetto
import h5py
import legenddataflowscripts as ldfs
import lgdo
import numpy as np
from dbetto import AttrsDict, TextDB
from legendmeta import LegendMetadata
from numpy.typing import ArrayLike
from reboost.hpge.psd import _current_pulse_model as current_pulse_model

from . import SimflowConfig, nersc

log = logging.getLogger(__name__)


def _merge_defaults(user: dict, default: dict) -> dict:
    """Recursively merge default values into user configuration.

    Merges values from `default` into `user` without overwriting existing
    user values. For nested dictionaries, performs recursive merge.

    Parameters
    ----------
    user
        User configuration dictionary.
    default
        Default configuration dictionary.

    Returns
    -------
    dict
        Merged configuration dictionary with user values taking precedence.
    """
    result = dict(default)
    for k, v in user.items():
        if k in result and isinstance(result[k], dict) and isinstance(v, dict):
            result[k] = _merge_defaults(v, result[k])
        else:
            result[k] = v
    return result


def _make_path(d):
    for k, v in d.items():
        if isinstance(v, str):
            d[k] = Path(v).resolve()
        else:
            d[k] = _make_path(v)
    return d


def init_simflow_context(raw_config: dict, workflow=None) -> AttrsDict:
    """Pre-process and sanitize the Simflow configuration.

    - set default configuration fields;
    - substitute ``$_`` and environment variables;
    - convert to :class:`~dbetto.attrsdict.AttrsDict`;
    - cast filesystem paths to :class:`pathlib.Path`;
    - clone and configure `legend-metadata`;
    - attach a :class:`legendmeta.LegendMetadata` instance to the Simflow
      configuration;
    - export important environment variables.

    Parameters
    ----------
    raw_config
        path to the Simflow configuration file.
    workflow
        Snakemake workflow instance. If None, occurrences of ``$_`` in the
        configuration will be replaced with the path to the current working
        directory.

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

    if workflow is None:
        ldfs.subst_vars(
            raw_config,
            var_values={"_": Path().resolve()},
            use_env=True,
            ignore_missing=False,
        )
    else:
        ldfs.workflow.utils.subst_vars_in_snakemake_config(workflow, raw_config)

    config = AttrsDict(raw_config)

    # convert all strings in the "paths" block to pathlib.Path
    config["paths"] = _make_path(config.paths)

    if "l200data" in config.paths:
        config["paths"]["l200data"] = nersc.dvs_ro(config, config.paths.l200data)

    # NOTE: this will attempt a clone of legend-metadata, if the directory does not exist
    metadata = LegendMetadata(config.paths.metadata, lazy=True)

    if "legend_metadata_version" in config:
        metadata.checkout(config.legend_metadata_version)

    # NOTE: read only path on NERSC, we are not going to modify the db
    # NOTE: don't use lazy=True, we need a fully functional TextDB
    config["metadata"] = LegendMetadata(nersc.dvs_ro(config, config.paths.metadata))

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


def setup_logdir_link(config: SimflowConfig, proctime: str) -> None:
    """Set up the timestamp-tagged directory for the workflow log files.

    Parameters
    ----------
    config
        Simflow configuration object.
    proctime
        Processing time identifier for the log directory.
    """
    logdir = Path(config.paths.log)
    logdir.mkdir(parents=True, exist_ok=True)

    # create a handy link to access latest log directory
    link = logdir / "latest"
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(proctime, target_is_directory=True)


def lookup_dataflow_config(l200data: Path | str) -> AttrsDict:
    """Finds and loads the dataflow configuration file.

    Parameters
    ----------
    l200data
        The path to the L200 data production cycle.

    Returns
    -------
    the dataflow configuration file as a dictionary with substitutions
    performed.
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

    df_cfg = AttrsDict(dbetto.utils.load_dict(df_cfgs[0]))
    df_cfg = df_cfg.setups.l200 if ("setups" in df_cfg) else df_cfg

    # substitute vars
    ldfs.subst_vars(
        df_cfg,
        var_values={"_": l200data.resolve()},
        use_env=False,
        ignore_missing=True,
    )

    # convert all strings in the "paths" block to pathlib.Path
    df_cfg["paths"] = _make_path(df_cfg.paths)

    return df_cfg


def init_generated_pars_db(
    l200data: str | Path, tier: str | None = None, lazy: bool = True
) -> TextDB:
    """Initializes the pars database from a LEGEND-200 data production.

    Parameters
    ----------
    l200data
        path to LEGEND-200 data production cycle.
    tier
        pars subfolder referring to a `tier`. If None, return the full `par`
        database.
    lazy
        see :class:`~dbetto.textdb.TextDB`.
    """
    dataflow_config = lookup_dataflow_config(l200data)
    return TextDB(
        dataflow_config.paths["par" if tier is None else f"par_{tier}"], lazy=lazy
    )


def get_hit_tier_name(l200data: str) -> str:
    """Extract the name of the hit tier for this production cycle.

    If the `pht` tier is present this is used else the `hit` tier is used.

    Parameters
    ----------
    l200data
        Path to the production cycle of l200 data.
    """

    df_cfg = lookup_dataflow_config(l200data).paths

    # first check if pht exists
    has_pht = ("tier_pht" in df_cfg) and Path(df_cfg.tier_pht).exists()
    has_hit = ("tier_hit" in df_cfg) and Path(df_cfg.tier_hit).exists()

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


def add_field_string(name: str, chunk: lgdo.Table, data: str) -> None:
    """Add a string to the output table.

    This is done in an HDF5-friendly way by storing the runid as a fixed-length
    string.
    """
    # NOTE: is this ok? will it compress well?
    dtype = h5py.string_dtype(encoding="utf-8", length=len(data))
    data_array = np.full(len(chunk), fill_value=data, dtype=dtype)
    chunk.add_field(name, lgdo.Array(data_array))


def check_nans_leq(array: ArrayLike, name: str, less_than_frac: float = 0.1) -> None:
    """Raise an exception if the fraction of NaN values in `array` is above threshold.

    Parameters
    ----------
    array
        the array to analyze.
    name
        array name for exception message.
    less_than_frac
        raise exception if fraction of NaNs is above this threshold.
    """
    flat = ak.ravel(array)
    n_el = len(flat)
    n_nans = ak.sum(ak.is_none(ak.nan_to_none(flat)))
    if (n_nans / n_el) > less_than_frac:
        msg = f"more than {100 * less_than_frac}% of NaNs detected in array {name}!"
        raise RuntimeError(msg)


def sorted_by(subset: Sequence, order: Sequence) -> list:
    """Sort a sequence according to a specified order, dropping duplicates."""
    pos = {k: i for i, k in enumerate(order)}

    # stable sort by order
    out = sorted(subset, key=pos.__getitem__)

    # deduplicate while preserving order
    seen = set()
    uniq = []
    for x in out:
        if x not in seen:
            seen.add(x)
            uniq.append(x)

    return uniq
