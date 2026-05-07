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
from numbers import Real
from pathlib import Path
from typing import Any

import awkward as ak
import dbetto
import h5py
import legenddataflowscripts as ldfs
import lgdo
import numpy as np
import yaml
from dbetto import AttrsDict, TextDB
from git.exc import GitCommandError
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


def apply_path_defaults(paths: dict) -> None:
    """Set default values for optional path keys derived from ``paths.pars``.

    The following keys are optional in the Simflow configuration and, if
    absent, are derived from ``paths.pars``:

    - ``geom``: defaults to ``{paths.pars}/geom``
    - ``dtmaps``: defaults to ``{paths.pars}/hpge/dtmaps``

    Parameters
    ----------
    paths
        The ``paths`` section of the Simflow configuration, with all values
        already converted to :class:`pathlib.Path` objects.

    """
    if "geom" not in paths:
        paths["geom"] = paths["pars"] / "geom"
    if "dtmaps" not in paths:
        paths["dtmaps"] = paths["pars"] / "hpge/dtmaps"


def link_external_paths(
    config: AttrsDict,
    workflow_basedir: str | Path,
    *,
    logger: logging.Logger | None = None,
) -> None:
    """Symlink user-overridden paths back into their default locations.

    When the user has manually overridden a ``paths.<key>`` entry in
    ``simflow-config.yaml`` to point outside the current production cycle
    (e.g. reusing the ``hit`` tier from another production), this function
    creates a symlink at the canonical default location pointing to the
    override. Snakemake rules keep reading ``config.paths.<key>``
    directly; the symlink only exists so the prod cycle's own
    ``generated/`` tree shows the external data in the standard layout.

    The default locations are computed from the simflow's own template at
    ``<workflow_basedir>/../templates/default.yaml``, with ``$_``
    substituted to the current working directory (the prod cycle root in
    the standard Snakemake invocation). Created symlinks are relative to
    the destination's parent, keeping the prod cycle portable.

    For each supported key:

    - if ``config.paths.<key>`` resolves to the default location, no
      override is in effect; any stale symlink at that location is removed;
    - otherwise a symlink is created (or refreshed) at the default
      location pointing to ``config.paths.<key>``.

    Real directories at the default location are never touched. The call
    is a safe no-op when ``<workflow_basedir>/../templates/default.yaml``
    does not exist.

    Supported keys (relative to ``paths``): every ``tier.<name>``,
    ``pars``, ``macros``, ``geom`` and ``dtmaps``. Default paths for
    ``geom`` and ``dtmaps`` fall back to ``<pars>/geom`` and
    ``<pars>/hpge/dtmaps`` when absent from the template (mirroring
    :func:`apply_path_defaults`).

    Parameters
    ----------
    config
        Simflow configuration as returned by :func:`init_simflow_context`.
    workflow_basedir
        Snakemake workflow basedir (``workflow.basedir`` in a Snakefile).
        Used only to locate the simflow's default template.
    logger
        Logger to use for status messages (e.g. the Snakemake logger when
        called from a Snakefile). Defaults to the module logger.

    """
    log_ = logger if logger is not None else log
    template = Path(workflow_basedir).parent / "templates" / "default.yaml"
    if not template.is_file():
        return

    default_cfg = yaml.safe_load(template.read_text())
    if not isinstance(default_cfg, dict):
        return
    ldfs.subst_vars(
        default_cfg,
        var_values={"_": str(Path.cwd().resolve())},
        ignore_missing=True,
    )
    default_paths_d = default_cfg.get("paths", {})
    # mirrors apply_path_defaults() — geom/dtmaps are commented out by default
    if "pars" in default_paths_d:
        default_paths_d.setdefault("geom", str(Path(default_paths_d["pars"]) / "geom"))
        default_paths_d.setdefault(
            "dtmaps", str(Path(default_paths_d["pars"]) / "hpge/dtmaps")
        )

    candidates: list[tuple[Any, Any]] = []
    for tier_name, current in config.paths.get("tier", {}).items():
        candidates.append((current, (default_paths_d.get("tier") or {}).get(tier_name)))
    for key in ("pars", "macros", "geom", "dtmaps"):
        if key in config.paths:
            candidates.append((config.paths[key], default_paths_d.get(key)))

    for current_str, default_str in candidates:
        if default_str is None:
            continue
        # normalize away the NERSC /dvs_ro read-only mount: the canonical
        # filesystem location is /global, and using /dvs_ro here would force
        # the relative symlink target up to / and back down through the mount.
        current = Path(nersc.undo_dvs_ro(str(current_str)))
        default = Path(nersc.undo_dvs_ro(str(default_str)))

        # compare absolute paths without resolving symlinks: otherwise a
        # symlink we just installed at `default` would make resolve() equal
        # current and look like a "no override".
        # `os.path.abspath` (not Path.resolve) on purpose: don't follow symlinks
        if os.path.abspath(current) == os.path.abspath(default):  # noqa: PTH100
            if default.is_symlink():
                default.unlink()
            continue

        # current path is an override: ensure a symlink at the default points to it
        if default.is_symlink():
            if (default.parent / default.readlink()).resolve() == current.resolve():
                continue
            default.unlink()
        elif default.exists():
            log_.warning(
                "%s is an existing non-symlink path; ignoring override -> %s",
                default,
                current,
            )
            continue

        if not current.exists():
            log_.warning(
                "override target %s does not exist; %s will be a broken symlink",
                current,
                default,
            )
        elif not current.is_dir():
            log_.warning(
                "override target %s is not a directory; linking anyway", current
            )

        rel = Path(os.path.relpath(current, start=default.parent))
        default.parent.mkdir(parents=True, exist_ok=True)
        default.symlink_to(rel, target_is_directory=True)


def init_simflow_context(
    raw_config: dict | AttrsDict | str | Path,
    workflow=None,
    logger: logging.Logger | None = None,
) -> AttrsDict:
    """Pre-process and sanitize the Simflow configuration.

    Returns a dictionary with useful objects to be used in the Simflow
    Snakefiles (i.e. the "context"):

    - set default configuration fields;
    - substitute ``$_`` and environment variables;
    - convert to :class:`~dbetto.attrsdict.AttrsDict`;
    - cast filesystem paths to :class:`pathlib.Path`;
    - clone and configure `legend-metadata`;
    - attach a :class:`~legendmeta.legendmetadata.LegendMetadata` instance to
      the Simflow configuration;
    - export important environment variables.

    Parameters
    ----------
    raw_config
        Simflow configuration mapping or path to a configuration file.
    workflow
        Snakemake workflow instance. If None, occurrences of ``$_`` in the
        configuration will be replaced with the path to the current working
        directory.
    logger
        Logger to use for status messages (e.g. the Snakemake logger when
        called from a Snakefile). Defaults to the module logger.

    """
    log_ = logger if logger is not None else log
    if not raw_config:
        msg = "you must set a config file with --configfile"
        raise RuntimeError(msg)

    basedir = Path(workflow.basedir) if workflow is not None else Path().resolve()

    if isinstance(raw_config, (str, Path)):
        if workflow is None:
            basedir = Path(raw_config).resolve().parent
        raw_config = dbetto.utils.load_dict(raw_config)

    if not isinstance(raw_config, (AttrsDict, dict)):
        msg = "raw_config must be a dict, AttrsDict, or path to a config file"
        raise TypeError(msg)

    if isinstance(raw_config, AttrsDict) and isinstance(
        raw_config.get("metadata"), LegendMetadata
    ):
        config = raw_config
    else:
        raw_config = _merge_defaults(
            dict(raw_config),
            {
                "benchmark": {"enabled": False},
                "nersc": {"dvs_ro": False, "scratch": False},
            },
        )

        if workflow is None:
            ldfs.subst_vars(
                raw_config,
                var_values={"_": basedir},
                use_env=True,
                ignore_missing=False,
            )
        else:
            ldfs.workflow.utils.subst_vars_in_snakemake_config(workflow, raw_config)

        config = AttrsDict(raw_config)

        # convert all strings in the "paths" block to pathlib.Path
        config["paths"] = _make_path(config.paths)

        # fill in optional paths derived from paths.pars
        apply_path_defaults(config["paths"])

        if "l200data" in config.paths:
            config["paths"]["l200data"] = nersc.dvs_ro(config, config.paths.l200data)

        # NOTE: this will attempt a clone of legend-metadata, if the directory does not exist
        metadata = LegendMetadata(config.paths.metadata, lazy=True)

        if "legend_metadata_version" in config:
            log_.info(
                "checking out legend-metadata version %s",
                config.legend_metadata_version,
            )
            try:
                metadata.checkout(config.legend_metadata_version)
            except GitCommandError as e:
                log_.warning("could not checkout legend-metadata version: %s", e)

        # NOTE: read only path on NERSC, we are not going to modify the db
        # NOTE: don't use lazy=True, we need a fully functional TextDB
        config["metadata"] = LegendMetadata(nersc.dvs_ro(config, config.paths.metadata))

    # make sure all simflow plots are made with a consistent style
    # I have verified only that this variable is visible in scripts (not shell directives)
    mpl_config_dir = (
        Path(__file__).resolve().parent
        if workflow is None
        else basedir / "src/legendsimflow"
    )
    os.environ["MPLCONFIGDIR"] = str(mpl_config_dir)

    if "_proctime" not in config:
        benchmark_enabled = bool(
            config.get("benchmark", AttrsDict({"enabled": False})).get("enabled", False)
        )
        proctime = (
            "benchmark"
            if benchmark_enabled
            else datetime.now().strftime("%Y%m%dT%H%M%SZ")
        )
        config["_proctime"] = proctime

    return AttrsDict({"config": config, "basedir": basedir})


def setup_logdir_link(config: SimflowConfig) -> None:
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
    link.symlink_to(config._proctime, target_is_directory=True)


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


def get_evt_tier_name(l200data: str) -> str:
    """Extract the name of the evt tier for this production cycle.

    If the `pet` tier is present this is used else the `evt` tier is used.

    Parameters
    ----------
    l200data
        Path to the production cycle of l200 data.

    """
    df_cfg = lookup_dataflow_config(l200data).paths

    # first check if pet exists
    has_pet = ("tier_pet" in df_cfg) and Path(df_cfg.tier_pet).exists()
    has_evt = ("tier_evt" in df_cfg) and Path(df_cfg.tier_evt).exists()

    if has_pet:
        return "pet"
    if has_evt:
        return "evt"

    msg = f"The l200data {l200data} does not contain a valid pet or evt tier"
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
    """Get the :func:`scipy.optimize.curve_fit` parameter results as a dictionary."""
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


def check_nans_leq(
    array: ArrayLike,
    name: str,
    less_than_frac: float = 0.1,
    min_entries: int = 100,
) -> None:
    """Raise an exception if the fraction of NaN values in `array` is above threshold.

    Parameters
    ----------
    array
        the array to analyze.
    name
        array name for exception message.
    less_than_frac
        raise exception if fraction of NaNs is above this threshold.
    min_entries
        minimum number of entries required to apply the fraction check. With
        fewer entries, a warning is logged instead of raising an exception.

    """
    flat = ak.ravel(array)
    n_el = len(flat)
    n_nans = ak.sum(ak.is_none(ak.nan_to_none(flat)))
    if n_el < min_entries:
        if n_nans > 0:
            msg = f"{n_nans}/{n_el} NaNs in {name}, but too few entries to apply fraction check"
            log.warning(msg)
        return
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


def sanitize_dict_with_defaults(read_dict: dict, defaults: dict) -> dict:
    """Swap-in defaults when values are illegal."""
    out = read_dict.copy()

    for key, sub in defaults.items():
        out.setdefault(key, {})
        for field, default_val in sub.items():
            val = out[key].get(field, default_val)

            ok = isinstance(val, type(default_val))
            if isinstance(default_val, Real):
                ok = isinstance(val, Real)

            if not ok:
                msg = (
                    f"{key}.{field}={val!r} has wrong type; "
                    f"using default {default_val!r}"
                )
                log.warning(msg)
                out[key][field] = default_val
            else:
                out[key][field] = val

    return out


def get_dict_value(d: dict, field: str, default: Any | None = None) -> Any:
    """Return a value from a nested dictionary using a dot-separated field path.

    Parameters
    ----------
    d
        Dictionary to query.
    field
        Dot-separated path (e.g. ``"a.b.c"``).
    default
        Value returned if the field is not found. Defaults to ``None``.

    """
    _ptr = d
    try:
        for segment in field.split("."):
            _ptr = _ptr[segment]
    except (KeyError, TypeError):
        return default

    return _ptr
