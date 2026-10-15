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

import h5py
from dbetto import AttrsDict, TextDB
from dbetto.catalog import Catalog
from legendmeta import LegendMetadata

from . import SimflowConfig, aggregate, patterns, utils
from .exceptions import SimflowConfigError
from .metadata import (
    get_par_settings,
    get_runlist,
    get_tier_settings,
    parse_runid,
    reference_cal_run,
    runinfo,
    simpars,
    simulates_hpge_psd,
)

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

#: Objects that must exist in the LH5 files of the data production.
#: ``{rawid}`` is replaced with the ``daq.rawid`` of an HPGe detector. The keys
#: are:
#:
#: - ``hpge_hit``: hit-tier table read to fit the current-pulse model;
#: - ``hpge_raw``: raw-tier waveforms read for the current-pulse model and the
#:   superpulses;
#: - ``evt_superpulses``: event-tier fields read to select the superpulse
#:   events. The time fields set in the ``superpulses`` settings are added;
#: - ``evt_rc_forced_trigger``, ``evt_rc_noise_trigger``: event-tier fields
#:   read to add random coincidences, for each ``random_coincidence_mode``;
#: - ``evt_impurity_tuning``: event-tier fields read to tune the HPGe
#:   impurities.
#:
#: The format is described in :func:`check_spec`.
LH5_SPECS = {
    "hpge_hit": {
        "ch{rawid}/hit": {
            "cuspEmax_ctc_cal": {},
            "cuspEmax_cal": {},
            "AoE_Classifier": {},
            "dt_eff": {},
        }
    },
    "hpge_raw": {
        "ch{rawid}/raw": {
            "waveform_windowed": {"values": {}, "t0": {}, "dt": {}},
            "waveform_presummed": {"values": {}, "t0": {}, "dt": {}},
        }
    },
    "evt_superpulses": {
        "evt": {
            "coincident/puls": {},
            "geds/energy": {},
            "geds/hit_idx": {},
            "geds/multiplicity": {},
            "geds/rawid": {},
            "geds/quality/is_good_channel": {},
            "geds/quality/is_bb_like": {},
            "geds/psd/low_aoe/value": {},
            "geds/psd/high_aoe/value": {},
        }
    },
    "evt_rc_forced_trigger": {
        "evt": {
            "trigger/is_forced": {},
            "coincident/geds": {},
            "coincident/muon_offline": {},
            "coincident/puls": {},
            "spms/energy": {},
            "spms/t0": {},
            "spms/rawid": {},
        }
    },
    "evt_rc_noise_trigger": {
        "evt": {
            "spms/energy": {},
            "spms/t0": {},
            "spms/rawid": {},
            "spms/is_trig_coin_pulse": {},
            "spms/quality/is_physical": {},
        }
    },
    "evt_impurity_tuning": {
        "evt": {
            "trigger/is_forced": {},
            "coincident/puls": {},
            "coincident/muon": {},
            "coincident/muon_offline": {},
            "geds/energy": {},
            "geds/detector_name": {},
            "geds/psd/low_aoe/time": {},
            "geds/quality/is_good_channel": {},
            "geds/quality/is_bb_like": {},
            "spms/event_t0": {},
            "spms/energy_sum": {},
        }
    },
}


def check_spec(
    obj: Mapping | TextDB, spec: Mapping[str, Mapping | type[Path]], where: str = ""
) -> list[str]:
    """Compare a dictionary or a text database with a spec.

    The spec is a dictionary. Each key must be found in `obj`, and its value is
    compared with the spec value:

    - ``{}``: the key must exist, its value can be anything;
    - a non-empty dictionary: the value must be a dictionary (or a folder of a
      text database, or a group of an HDF5 file) and is compared with it
      recursively;
    - :class:`~pathlib.Path`: the value must be a path to an existing file or
      folder.

    A key can be a path through the nested dictionaries or folders, like
    ``datasets/runinfo``. Alternatives are separated by ``|``, like
    ``tier_pht|tier_hit``: at least one of them must match the spec.

    Use a lazy :class:`~dbetto.textdb.TextDB` to check files and folders: only
    the files named in the spec are read. An open :class:`h5py.File` can be
    checked too: only the names of its objects are read, not the data.

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
                value = obj
                for part in alt.split("/"):
                    value = value[part]
            except (KeyError, TypeError, ValueError, FileNotFoundError):
                errs = [f"{loc}: not found"]
            else:
                if subspec is Path:
                    exists = isinstance(value, Path) and value.exists()
                    errs = [] if exists else [f"{loc}: {value} does not exist"]
                elif not subspec:
                    errs = []
                elif isinstance(value, Mapping | TextDB):
                    errs = check_spec(value, subspec, loc)
                else:
                    errs = [f"{loc}: must be a dictionary or a folder"]

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


def _skip(config: SimflowConfig) -> bool:
    """Whether ``paths.l200data`` is not set or ``validate_l200data`` is false."""
    if "l200data" not in config.paths:
        return True
    if not config.get("validate_l200data", True):
        log.warning(
            "validate_l200data is false: the data production at %s is not validated",
            config.paths.l200data,
        )
        return True
    return False


def _data_problems(
    config: SimflowConfig, l200data: Path, df_cfg: AttrsDict
) -> list[str]:
    """Check the parameter and LH5 files that the Simflow reads from the data production."""
    meta = config.metadata
    hit_settings = get_tier_settings(config, "hit")
    evt_settings = get_tier_settings(config, "evt")
    runids = sorted(aggregate.gen_list_of_all_runids(config))
    hit_tier = utils.get_hit_tier_name(l200data)
    evt_tier = utils.get_evt_tier_name(l200data)
    problems = []

    # the parameter files are large: find them without reading them
    pars_dir = df_cfg.paths[f"par_{hit_tier}"]
    validity = [
        f
        for ext in (".yaml", ".yml", ".json", ".jsonl")
        if (f := pars_dir / f"validity{ext}").is_file()
    ]
    if len(validity) != 1:
        problems.append(f"{pars_dir}: expected one validity file")
    else:
        catalog = Catalog.read_from(validity[0])
        for runid in runids:
            tstamp = runinfo(meta, runid).start_key
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

    # runs to check for each (tier, LH5_SPECS key)
    checks: dict[tuple[str, str], set[str]] = {}
    specs = dict(LH5_SPECS)

    tune_impurities = hit_settings.get("tune_hpge_impurities_on_data", False)
    currmod = simpars(meta, "geds.currmod", runids[0], config.experiment, default={})
    currmod_from_data = simulates_hpge_psd(config) and currmod.get("default") is None
    superpulses_from_data = tune_impurities or (
        hit_settings.get("simulate_psd_with_psl", False)
        and patterns.compute_superpulses(config, runid=runids[0])
    )
    if currmod_from_data or superpulses_from_data:
        cal_runs = {r: reference_cal_run(meta, r) for r in runids}

    if currmod_from_data:
        runs = set(runids) | set(cal_runs.values())
        checks.setdefault((hit_tier, "hpge_hit"), set()).update(runs)
        checks.setdefault(("raw", "hpge_raw"), set()).update(runs)

    if superpulses_from_data:
        # heavy imports (lh5, matplotlib), needed only here
        from .superpulses import DEFAULT_SETTINGS  # noqa: PLC0415

        settings = DEFAULT_SETTINGS | get_par_settings(config, "superpulses")
        tier = evt_tier if settings["tier_name"] == "auto" else settings["tier_name"]
        runs = {cal_runs[r] if r.endswith("-phy") else r for r in runids}
        if tier in ("evt", "pet"):
            time_fields = [settings["t0_field"], settings["end_time_field"]]
            specs["evt_superpulses"] = {
                "evt": LH5_SPECS["evt_superpulses"]["evt"]
                | {f: {} for f in time_fields if f}
            }
            checks.setdefault((tier, "evt_superpulses"), set()).update(runs)
        else:
            checks.setdefault((tier, ""), set()).update(runs)
        checks.setdefault(("raw", "hpge_raw"), set()).update(runs)

    if evt_settings.get("add_random_coincidences", False):
        rc_runid = evt_settings.get("random_coincidence_runid")
        if isinstance(rc_runid, Mapping):
            runs = set(rc_runid.values())
        else:
            runs = {rc_runid} if rc_runid is not None else set(runids)
        mode = evt_settings.get("random_coincidence_mode", "forced_trigger")
        checks.setdefault((evt_tier, f"evt_rc_{mode}"), set()).update(runs)

    if tune_impurities:
        runs = {
            r
            for simid in aggregate.gen_list_of_tuning_simids(config)
            for r in get_runlist(config, simid)
        }
        checks.setdefault((evt_tier, "evt_impurity_tuning"), set()).update(runs)

    if not checks:
        return problems

    # the data production formats all runs and detectors in the same way: check
    # that each run has files, but the content of only one file and detector
    chmap = meta.channelmap(runinfo(meta, runids[0]).start_key)
    rawid = next(
        (
            ch.daq.rawid
            for ch in chmap.values()
            if ch.get("system") == "geds"
            and ch.get("analysis", {}).get("usability") == "on"
        ),
        None,
    )
    if rawid is None:
        return [
            *problems,
            f"{runids[0]}: no HPGe detector is on in the Simflow metadata",
        ]

    first_files: dict[Path, Path | None] = {}
    for (tier, name), runs in checks.items():
        first = None
        for runid in sorted(runs):
            _, period, run, datatype = parse_runid(runid)
            folder = (
                df_cfg.paths[f"tier_{tier}"]
                / datatype
                / f"p{period:02d}"
                / f"r{run:03d}"
            )
            if folder not in first_files:
                first_files[folder] = next(folder.glob("*.lh5"), None)
                if first_files[folder] is None:
                    problems.append(f"{runid}: no LH5 files in {folder}")
            first = first or first_files[folder]

        if first is not None and name:
            with h5py.File(first, "r") as f:
                problems += check_spec(
                    f,
                    {k.format(rawid=rawid): v for k, v in specs[name].items()},
                    str(first),
                )

    return problems


def l200data_layout(config: SimflowConfig) -> None:
    """Check that the data production at ``paths.l200data`` has the expected format.

    Checks:

    - :data:`L200DATA_SPEC`: the metadata in ``inputs/``;
    - :data:`DATAFLOW_CONFIG_SPEC`: the dataflow configuration file and the
      tier and parameter folders it points to;
    - a hit-tier parameter file is valid at the start of each run of the
      Simflow;
    - the ICPC DSP configuration holds the current smoothing settings (see
      :func:`.utils.lookup_curr_mw_pars`);
    - the LH5 files required by the Simflow settings exist for each run, and
      their first file has the objects in :data:`LH5_SPECS`. Only one HPGe
      detector is checked: the first one that is on in the first run of the
      Simflow.

    All problems are reported at once. Does nothing if ``paths.l200data`` is
    not set, or if ``validate_l200data`` is false.

    Parameters
    ----------
    config
        Simflow configuration object.
    """
    if _skip(config):
        return

    l200data = Path(config.paths.l200data)
    if not l200data.is_dir():
        msg = f"the data production {l200data} is not a folder"
        raise SimflowConfigError(msg)

    problems = check_spec(TextDB(l200data, lazy=True), L200DATA_SPEC, str(l200data))

    try:
        df_cfg = utils.lookup_dataflow_config(l200data)
        df_problems = check_spec(df_cfg, DATAFLOW_CONFIG_SPEC, "dataflow config")
    except RuntimeError as e:
        df_problems = [str(e)]
    problems += df_problems

    # the files can be checked only if their folders are known
    if not df_problems:
        problems += _data_problems(config, l200data, df_cfg)

    try:
        utils.lookup_curr_mw_pars(l200data)
    except (RuntimeError, KeyError, ValueError) as e:
        problems.append(str(e))

    if problems:
        msg = (
            f"the data production {l200data} does not have the expected format:\n"
            + "\n".join(f"  {p}" for p in problems)
        )
        raise SimflowConfigError(msg)


def detinfo_against_l200data(
    config: SimflowConfig,
    detinfo: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> None:
    """Check that the detector usabilities match the data production metadata.

    Does nothing if ``paths.l200data`` is not set, or if ``validate_l200data``
    is false.

    Parameters
    ----------
    config
        Simflow configuration object.
    detinfo
        Detector info computed from the Simflow metadata, as returned by
        :func:`.aggregate.pivot_detinfo`.
    """
    if _skip(config):
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
        raise SimflowConfigError(msg)
