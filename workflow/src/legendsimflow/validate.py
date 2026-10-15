# Copyright (C) 2026 Luigi Pertoldi <gipert@pm.me>
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

"""Checks of the Simflow configuration against the LEGEND-200 data production.

The data production is found at ``paths.l200data`` in the Simflow
configuration. Each check does nothing if this path is not set.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path

from dbetto import AttrsDict, TextDB
from dbetto.catalog import Catalog
from legendmeta import LegendMetadata

from . import SimflowConfig, aggregate, utils
from .exceptions import SimflowConfigError
from .metadata import runinfo

log = logging.getLogger(__name__)

#: Files and folders that must exist in ``paths.l200data``. The format is
#: described in :func:`check_spec`.
L200DATA_SPEC = {
    "inputs": {
        "datasets/runinfo|dataprod/runinfo": {},
        "datasets/statuses|dataprod/config": {},
        "hardware/configuration/channelmaps": {},
        "hardware/detectors": {},
    }
}

#: Required content of the dataflow configuration in ``paths.l200data``, as
#: returned by :func:`.utils.lookup_dataflow_config`. The format is described
#: in :func:`check_spec`.
DATAFLOW_CONFIG_SPEC = {
    "paths": {
        "tier_raw": Path,
        "tier_pht|tier_hit": Path,
        "tier_pet|tier_evt": Path,
        "par_pht|par_hit": Path,
    }
}


def check_spec(
    obj: Mapping | TextDB, spec: Mapping[str, Mapping | type], where: str = ""
) -> list[str]:
    """Compare a dictionary or a text database with a spec.

    The spec is a dictionary. Each key must be found in `obj`, and its value is
    compared with the spec value:

    - a dictionary: the value must be a dictionary (or a folder of a text
      database) and is compared with it recursively. Use ``{}`` to only
      require that the key exists;
    - :class:`~pathlib.Path`: the value must be a path to an existing file or
      folder;
    - another type, e.g. ``str``: the value must be of this type (``float``
      also accepts integers).

    A key can be a path through the nested dictionaries or folders, like
    ``datasets/runinfo``. Alternatives are separated by ``|``, like
    ``tier_pht|tier_hit``: at least one of them must match the spec.

    Use a lazy :class:`~dbetto.textdb.TextDB` to check files and folders: only
    the files named in the spec are read.

    Parameters
    ----------
    obj
        dictionary or text database to check.
    spec
        the spec.
    where
        location of `obj`, used as prefix in the messages.

    Returns
    -------
    one message for each problem found. The list is empty if `obj` matches
    the spec.
    """
    problems = []
    for key, subspec in spec.items():
        alternatives = []
        for alt in key.split("|"):
            loc = f"{where}/{alt}" if where else alt
            try:
                if isinstance(obj, TextDB):
                    value = obj[alt]
                else:
                    value = obj
                    for part in alt.split("/"):
                        value = value[part]
            except (KeyError, TypeError, FileNotFoundError):
                alternatives.append([f"{loc}: not found"])
                continue

            if isinstance(subspec, Mapping):
                if isinstance(value, Mapping | TextDB):
                    errs = check_spec(value, subspec, loc)
                else:
                    errs = [f"{loc}: must be a dictionary or a folder"]
            elif subspec is Path:
                exists = isinstance(value, Path) and value.exists()
                errs = [] if exists else [f"{loc}: {value} does not exist"]
            elif not isinstance(value, subspec) and not (
                subspec is float and isinstance(value, int)
            ):
                errs = [
                    f"{loc}: must be of type {subspec.__name__}, "
                    f"found {type(value).__name__}"
                ]
            else:
                errs = []

            if not errs:
                break
            alternatives.append(errs)
        else:
            if len(alternatives) > 1:
                problems.append(f"no valid alternative among {key!r} in {where!r}:")
                problems += [f"  {err}" for errs in alternatives for err in errs]
            else:
                problems += alternatives[0]

    return problems


def _report(config: SimflowConfig, msg: str) -> None:
    """Raise an error, or log a warning if ``allow_invalid_l200data`` is true."""
    if config.get("allow_invalid_l200data", False):
        log.warning(msg)
    else:
        raise SimflowConfigError(msg)


def l200data_layout(config: SimflowConfig) -> None:
    """Check that the data production at ``paths.l200data`` has the expected format.

    Checks:

    - :data:`L200DATA_SPEC`: the metadata in ``inputs/``;
    - :data:`DATAFLOW_CONFIG_SPEC`: the dataflow configuration file and the
      tier and parameter folders it points to;
    - a hit-tier parameter file is valid at the start of each run of the
      Simflow;
    - the ICPC DSP configuration holds the current smoothing settings (see
      :func:`.utils.lookup_curr_mw_pars`).

    All problems are reported at once. Does nothing if ``paths.l200data`` is
    not set. Logs a warning instead of raising an error if
    ``allow_invalid_l200data`` is true.

    Parameters
    ----------
    config
        Simflow configuration object.
    """
    if "l200data" not in config.paths:
        return

    l200data = Path(config.paths.l200data)
    if not l200data.is_dir():
        _report(config, f"the data production {l200data} is not a folder")
        return

    problems = check_spec(TextDB(l200data, lazy=True), L200DATA_SPEC, str(l200data))

    try:
        df_cfg = utils.lookup_dataflow_config(l200data)
    except RuntimeError as e:
        problems.append(str(e))
    else:
        df_problems = check_spec(df_cfg, DATAFLOW_CONFIG_SPEC, "dataflow config")
        problems += df_problems

        if not df_problems:
            # the parameter files are large: find them without reading them
            pars_dir = df_cfg.paths[f"par_{utils.get_hit_tier_name(l200data)}"]
            validity = list(pars_dir.glob("validity.*"))
            if len(validity) != 1:
                problems.append(f"{pars_dir}: expected one validity file")
            else:
                catalog = Catalog.read_from(validity[0])
                for runid in sorted(aggregate.gen_list_of_all_runids(config)):
                    tstamp = runinfo(config.metadata, runid).start_key
                    try:
                        files = catalog.valid_for(tstamp)
                    except RuntimeError:
                        files = []
                    if not files:
                        problems.append(
                            f"{runid}: no parameter file in {pars_dir} is valid at {tstamp}"
                        )
                    problems += [
                        f"{runid}: {pars_dir / f} does not exist"
                        for f in files
                        if not (pars_dir / f).is_file()
                    ]

    try:
        utils.lookup_curr_mw_pars(l200data)
    except (RuntimeError, KeyError, ValueError) as e:
        problems.append(str(e))

    if problems:
        _report(
            config,
            f"the data production {l200data} does not have the expected format:\n"
            + "\n".join(f"  {p}" for p in problems),
        )


def detinfo_against_l200data(
    config: SimflowConfig,
    detinfo: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> None:
    """Check that the detector usabilities match the data production metadata.

    Does nothing if ``paths.l200data`` is not set. Logs a warning instead of
    raising an error if ``allow_invalid_l200data`` is true.

    Parameters
    ----------
    config
        Simflow configuration object.
    detinfo
        Detector info computed from the Simflow metadata, as returned by
        :func:`.aggregate.pivot_detinfo`.
    """
    if "l200data" not in config.paths:
        return

    inputs = Path(config.paths.l200data) / "inputs"
    data_config = dict(config)
    data_config["metadata"] = LegendMetadata(str(inputs), lazy=True)
    data_config = AttrsDict(data_config)

    flags = ["usability", "psd_usability"]
    # the data production metadata has no simprod config: query the same runs
    runids = {runid for flag in flags for runid in detinfo.get(flag, {})}
    data_detinfo = aggregate.pivot_detinfo(
        aggregate.gen_list_of_all_usabilities(data_config, runids).to_dict()
    )

    for flag in flags:
        sim = detinfo.get(flag, {})
        data = data_detinfo.get(flag, {})
        if sim == data:
            continue

        diffs = []
        for runid in sorted(set(sim) | set(data)):
            sim_run = sim.get(runid, {})
            data_run = data.get(runid, {})
            for det in sorted(set(sim_run) | set(data_run)):
                sim_val = sim_run.get(det, "<missing>")
                data_val = data_run.get(det, "<missing>")
                if sim_val != data_val:
                    diffs.append(
                        f"  {runid} {det}: Simflow metadata has {sim_val!r}, "
                        f"data production has {data_val!r}"
                    )

        msg = (
            f"'{flag}' differs between the Simflow metadata and the data "
            f"production metadata in {inputs} for {len(diffs)} detector(s):\n"
            + "\n".join(diffs)
        )
        _report(config, msg)
