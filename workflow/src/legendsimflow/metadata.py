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

import logging
import re
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path

import lh5
from dbetto import AttrsDict
from dbetto.catalog import Catalog
from legendmeta import MetadataRepository
from legendmeta.police import validate_dict_schema
from snakemake.iocontainers import Wildcards

from . import SimflowConfig, defaults, utils
from .exceptions import SimflowConfigError

_MISSING = object()

log = logging.getLogger(__name__)

USABILITY_CODE = {
    "on": 0,
    "ac": 1,
    "off": 2,
}

PSD_USABILITY_CODE = {
    "valid": 0,
    "present": 1,
    "missing": 2,
}


def get_simconfig(
    config: SimflowConfig,
    tier: str,
    simid: str | None = None,
    field: str | None = None,
) -> AttrsDict:
    """Return the simulation configuration for the given tier and simid.

    Raise :class:`~legendsimflow.exceptions.SimflowConfigError` if any key is
    not found.

    Parameters
    ----------
    config
        Simflow configuration object.
    tier
        Tier name.
    simid
        Simulation identifier.
    field
        If not ``None``, return the value of this key in the simconfig.

    """
    try:
        _m = config.metadata.simprod.config
    except FileNotFoundError as e:
        raise SimflowConfigError(e) from e

    block = f"simprod.config.tier.{tier}.{config.experiment}.simconfig"
    try:
        if simid is None:
            block = f"simprod.config.tier.{tier}.{config.experiment}"
            simcfg = _m.tier[tier][config.experiment].simconfig
            validate_simconfig_keys(simcfg, block + ".simconfig")
            return simcfg
        if field is None:
            return _m.tier[tier][config.experiment].simconfig[simid]
        return _m.tier[tier][config.experiment].simconfig[simid][field]

    except KeyError as e:
        msg = f"key {e} not found!"
        raise SimflowConfigError(msg, block) from e
    except FileNotFoundError as e:
        raise SimflowConfigError(e, block) from e


def get_tier_settings(config: SimflowConfig, tier: str) -> AttrsDict:
    """Return the settings block for *tier* and the current experiment.

    Fields missing in the user settings (or the whole block, if absent) are
    taken from :data:`.defaults.TIER_SETTINGS`.
    """
    try:
        user = config.metadata.simprod.config.tier[tier][config.experiment].settings
    except FileNotFoundError:
        user = {}
    return defaults.merge(user, defaults.TIER_SETTINGS.get(tier, {}))


def simulates_hpge_psd(config: SimflowConfig) -> bool:
    """Whether HPGe pulse shapes are simulated: ``simulate_psd`` and no ``skip_hit``."""
    return (
        not get_tier_settings(config, "evt").skip_hit
        and get_tier_settings(config, "hit").simulate_psd
    )


def get_opt_scintillators(config: SimflowConfig) -> AttrsDict:
    """Return the light sources of the ``opt`` tier, ``<source name> -> settings``.

    Each source is a scintillator whose light reaches the SiPMs. Its settings
    are:

    - ``table``: name of its table in the ``stp`` files, i.e. the name of the
      registered scintillator volume, or the ntuple name shared by several
      volumes (e.g. ``pen`` for all PEN plates);
    - ``material``: scintillation model (yield and time profile), as known to
      reboost: ``lar``, ``pen`` or ``fiber``;
    - ``optmap_scaling_factor``: taken from the top-level value of the ``opt``
      tier settings, if not set for this source;
    - ``model``: ``None`` to use the optical map, else a mapping with the
      import path of a function that replaces the map, giving the detection
      probability of the photons emitted in each energy deposition
      (``function``), and the keyword arguments passed to it (``kwargs``,
      empty if not set).

    The optical map of source ``<name>`` is ``paths.optical_maps.<name>``.
    Without a ``scintillators`` block in the ``opt`` tier settings, the only
    source is ``lar``, with table ``scintillator_volume_name``.
    """
    settings = get_tier_settings(config, "opt")
    block = f"simprod.config.tier.opt.{config.experiment}.settings"

    sources = settings.scintillators
    if sources is None:
        if "scintillator_volume_name" not in settings:
            msg = "either 'scintillators' or 'scintillator_volume_name' must be set"
            raise SimflowConfigError(msg, block)
        sources = {
            "lar": {"table": settings.scintillator_volume_name, "material": "lar"}
        }
    elif not isinstance(sources, Mapping) or len(sources) == 0:
        msg = "'scintillators' must be a non-empty mapping <source name> -> settings"
        raise SimflowConfigError(msg, block)

    if "lar_patch" in sources:
        # paths.optical_maps.lar_patch is not the map of a light source
        msg = "'lar_patch' is not allowed as a name in 'scintillators'"
        raise SimflowConfigError(msg, block)

    sources = {name: dict(source) for name, source in sources.items()}
    tables = {}
    for name, source in sources.items():
        for key in ("table", "material"):
            if key not in source:
                msg = f"scintillators.{name} has no '{key}' field"
                raise SimflowConfigError(msg, block)
        if source["table"] in tables:
            msg = (
                f"scintillators.{name} and scintillators.{tables[source['table']]} "
                f"have the same table {source['table']}"
            )
            raise SimflowConfigError(msg, block)
        tables[source["table"]] = name

        model = source.setdefault("model", None)
        if model is not None:
            function = model.get("function") if isinstance(model, Mapping) else None
            if not isinstance(function, str) or "." not in function:
                msg = (
                    f"scintillators.{name}.model must be a mapping with the import "
                    "path of a function in 'function', e.g. 'mypackage.pen_light'"
                )
                raise SimflowConfigError(msg, block)
            if not isinstance(model.get("kwargs", {}), Mapping):
                msg = f"scintillators.{name}.model.kwargs must be a mapping"
                raise SimflowConfigError(msg, block)
            source["model"] = {"kwargs": {}} | dict(model)

        if "optmap_scaling_factor" in settings:
            source.setdefault("optmap_scaling_factor", settings.optmap_scaling_factor)

    return AttrsDict(sources)


def deferred_tier_setting(config: SimflowConfig, tier: str, key: str) -> Callable:
    """Make a Snakemake ``params:`` callable resolving a tier setting lazily."""
    return lambda wildcards: get_tier_settings(config, tier)[key]  # noqa: ARG005


def get_par_settings(
    config: SimflowConfig, par: str, default: Mapping | None = None
) -> AttrsDict:
    """Return the settings block for *par* and the current experiment.

    Fields missing in the user settings (or the whole block, if absent) are
    taken from `default`, or from :data:`.defaults.PAR_SETTINGS` if ``None``.
    """
    try:
        user = config.metadata.simprod.config.pars[config.experiment].geds[par].settings
    except FileNotFoundError:
        user = {}
    if default is None:
        default = defaults.PAR_SETTINGS.get(par, {})
    return defaults.merge(user, default)


def smk_hash_simconfig(
    config: SimflowConfig,
    wildcards: Wildcards,
    field: str | None = None,
    ignore: list | None = None,
    **kwargs,
) -> str:
    """Get the dictionary hash for use in Snakemake rules.

    Parameters
    ----------
    config
        Snakemake config.
    wildcards
        Snakemake wildcards object.
    field
        If not ``None``, return the value of this key in the simconfig.
    ignore
        Exclude these fields from the hash.
    kwargs
        provide a value for wildcards that might not be present in `wildcards`.

    """
    tier = kwargs["tier"] if "tier" in kwargs else wildcards.tier  # noqa: SIM401
    simid = kwargs["simid"] if "simid" in kwargs else wildcards.simid  # noqa: SIM401

    scfg = get_simconfig(config, tier, simid).copy()

    if field is not None:
        scfg = scfg.get(field)

    if ignore is not None:
        if not isinstance(ignore, tuple | list):
            ignore = [ignore]

        for f in ignore:
            if f in scfg:
                scfg.pop(f)

    return utils.hash_dict(scfg)


def extract_integer(file_path: Path) -> int:
    """Read a single integer from a file, stripping surrounding whitespace."""
    with file_path.open() as f:
        return int(f.read().strip())


def usability(
    metadata: MetadataRepository, det_name: str, runid: str, default: str | None = None
) -> str:
    """Get the usability for analysis of `det_name` in run `runid`.

    Looks for the ``analysis.usability`` metadata field in the channel map. By
    default, an error is thrown if no information is found. If `default` is set
    to a non-None value, it will be returned.
    """
    rinfo = runinfo(metadata, runid)
    chmap = metadata.channelmap(rinfo.start_key)
    if det_name in chmap and "analysis" in chmap[det_name]:
        return chmap[det_name].analysis.usability

    if default is None:
        msg = f"no usability metadata found for {det_name} in {runid} and no default provided"
        raise RuntimeError(msg)

    msg = (
        f"no usability metadata found for {det_name} in {runid}, returning {default!s}"
    )
    log.warning(msg)
    return default


def encode_usability(usability: str) -> int:
    """Encode the HPGe usability in an int."""
    return USABILITY_CODE[usability]


def decode_usability(usability_code: int) -> str:
    """Decode the HPGe usability (see :func:`encode_usability`)."""
    _codes = {v: k for k, v in USABILITY_CODE.items()}
    return _codes[usability_code]


def encode_psd_usability(psd_usability: str) -> int:
    """Encode the PSD usability in an int."""
    return PSD_USABILITY_CODE[psd_usability]


def decode_psd_usability(psd_usability_code: int) -> str:
    """Decode the PSD usability (see :func:`encode_psd_usability`)."""
    _codes = {v: k for k, v in PSD_USABILITY_CODE.items()}
    return _codes[psd_usability_code]


def parse_runid(runid: str) -> (str, int, int, str):
    """Extract `runid` fields.

    Returns the experiment, period, run and datatype as a tuple. Period and run
    are integers.
    """
    if not is_runid(runid):
        msg = f"{runid} is not a valid runid"
        raise ValueError(msg)

    experiment, period, run, datatype = re.split(r"\W+", runid)
    return experiment, int(period[1:]), int(run[1:]), datatype


def runinfo(metadata: MetadataRepository, runid: str) -> str:
    """Get the `datasets.runinfo` entry for a LEGEND run identifier.

    Parameters
    ----------
    metadata
        LEGEND metadata database.
    runid
        a run identifier in the format ``<experiment>-<period>-<run>-<datatype>``.

    """
    _, period, run, datatype = re.split(r"\W+", runid)
    try:
        runinfo_db = metadata.datasets.runinfo
    except FileNotFoundError:
        # old metadata releases kept it in dataprod
        runinfo_db = metadata.dataprod.runinfo
    return runinfo_db[period][run][datatype]


def reference_cal_run(metadata: MetadataRepository, runid: str) -> str:
    """The reference calibration run for `runid`.

    If an entry of the ``datasets.run_override`` validity file is valid at the
    start of `runid`, return the calibration run it lists (e.g.
    ``l200-p16-r006-cal-20251001T130640Z``). Otherwise return the calibration
    run of the same run or, if missing, of the closest previous run in the
    same period. A calibration `runid` is its own reference.
    """
    exp, period, run, datatype = re.split(r"\W+", runid)

    msg = f"looking for reference calibration run for {runid}"
    log.debug(msg)

    if datatype == "cal":
        return runid

    try:
        overrides = Catalog.build_catalog(metadata.datasets.run_override)
    except FileNotFoundError:
        overrides = None

    if overrides is not None:
        tstamp = runinfo(metadata, runid).start_key
        cal_keys = overrides.valid_for(tstamp, allow_none=True) or []
        if len(cal_keys) > 1:
            msg = f"more than one calibration run override is valid for {runid}"
            raise RuntimeError(msg)
        if cal_keys:
            cal_runid = cal_keys[0].rsplit("-", 1)[0]
            msg = f"calibration run override found for {runid}: {cal_runid}"
            log.debug(msg)
            return cal_runid

    p_runinfo = metadata.datasets.runinfo[period]

    if "cal" in p_runinfo[run]:
        return f"{exp}-{period}-{run}-cal"

    runs = sorted(p_runinfo.keys())
    index = runs.index(run)

    while True:
        index -= 1
        if index < 0:
            msg = f"there is no previous calibration run in {period} for {runid}"
            raise RuntimeError(msg)

        prev_r = runs[index]
        if "cal" in p_runinfo[prev_r]:
            return f"{exp}-{period}-{prev_r}-cal"


def simpars(
    metadata: MetadataRepository,
    par: str,
    runid: str,
    experiment: str,
    default: object = _MISSING,
) -> AttrsDict:
    """Extract simflow parameters for a certain LEGEND run.

    Queries the simflow parameters database stored under
    ``simprod.config.pars`` by experiment name `experiment`, parameter name
    `par` and LEGEND run identifier `runid`.

    Parameters
    ----------
    metadata
        LEGEND metadata database.
    par
        name of directory under ``metadata.simprod.config.pars.{experiment}``.
        Can be a nested property, as in e.g. ``geds.eresmod``. ``.`` and
        ``/`` are allowed separators.
    runid
        a run identifier in the format ``<experiment>-<period>-<run>-<datatype>``.
    experiment
        experiment identifier (e.g. ``l200cfg01``, ``l1000dsg01``). Selects
        the experiment-level subdirectory under ``simprod/config/pars/``.
    default
        value to return when the parameter directory is not found in the
        database or no validity entry matches `runid`. If not provided, such
        cases raise ``KeyError`` or ``LookupError``. Other errors (e.g.
        malformed YAML) are always re-raised regardless of this argument.

    """
    par = par.replace(".", "/")
    datatype = re.split(r"\W+", runid)[-1]
    try:
        directory = metadata["simprod/config/pars"][experiment][par]
        # raises FileNotFoundError if the folder holds only a settings.yaml
        directory["validity"]
        return directory.on(runinfo(metadata, runid).start_key, category=datatype)
    except (KeyError, LookupError, FileNotFoundError):
        if default is _MISSING:
            raise
        return default


def get_vtx_simconfig(config: SimflowConfig, simid: str) -> AttrsDict:
    """Get the vertex generation configuration for a stp-tier `simid`.

    Returns the ``vtx``-tier generator requested by the ``stp``-tier simulation
    with identifier `simid`.

    Parameters
    ----------
    config
        Snakemake config.
    simid
        simulation identifier.

    """
    vtx_key = set()
    sconfig = get_simconfig(config, "stp", simid)
    for field in ("generator", "confinement"):
        if field in sconfig and sconfig[field].startswith("~vertices:"):
            vtx_key.add(sconfig[field].partition(":")[2])

    if len(vtx_key) == 0:
        msg = f"{simid} does not specify vertices? This is unexpected, please file a bug report"
        raise RuntimeError(msg)

    if len(vtx_key) > 1:
        raise NotImplementedError()

    return get_simconfig(config, "vtx", vtx_key.pop())


def get_sanitized_fccd(metadata: MetadataRepository, det_name: str) -> float:
    """Return the FCCD value for `det_name`, falling back to 1 mm if the FCCD field is absent.

    Parameters
    ----------
    metadata
        LEGEND metadata database.
    det_name
        Detector name.

    """
    try:
        det_meta = metadata.hardware.detectors.germanium.diodes[det_name]
    except (FileNotFoundError, KeyError):
        msg = f"{det_name} diode metadata not found, setting FCCD to 1 mm"
        log.warning(msg)
        return 1.0

    has_fccd_meta = validate_dict_schema(
        det_meta.characterization,
        {"combined_0vbb_analysis": {"fccd_in_mm": {"value": 0.0}}},
        greedy=False,
        verbose=False,
    )

    if not has_fccd_meta:
        msg = f"{det_name} metadata does not seem to contain usable FCCD data, setting to 1 mm"
        log.warning(msg)
        fccd = 1
    else:
        fccd = det_meta.characterization.combined_0vbb_analysis.fccd_in_mm.value

    return fccd


def is_runid(runid: str) -> bool:
    """Whether a runid (run identifier) is correctly formatted.

    It should be in the form
    ``<experiment>-<period>-<run>-<datatype>``/``XXX-pNN-rMMM-AAA``
    where ``XXX`` is any alphanumeric experiment identifier.
    """
    return re.match(r"^[A-Za-z0-9]+-p(\d{2})-r(\d{3})-([A-Za-z]+)$", runid) is not None


def is_simid(simid: str) -> bool:
    r"""Whether a simid (simulation identifier) is correctly formatted.

    A valid simid must consist entirely of word characters (letters, digits,
    underscores) and hyphens, matching the pattern ``[-\w]+``. Dots and other
    special characters are not allowed; in particular, dots are forbidden
    because they are used as the delimiter in the simlist format
    ``<tier>.<simid>``.
    """
    return re.fullmatch(r"[-\w]+", simid, flags=re.ASCII) is not None


def validate_simconfig_keys(simconfig: Mapping, block: str | None = None) -> None:
    """Validate that all top-level keys of `simconfig` are valid simids.

    Raises :class:`~legendsimflow.exceptions.SimflowConfigError` listing every
    invalid key if any are found.

    Parameters
    ----------
    simconfig
        Dictionary whose top-level keys are expected to be simids (as loaded
        from a ``simconfig.yaml`` file).
    block
        Optional config block label included in the error message for context.
    """
    invalid_keys = [k for k in simconfig if not (isinstance(k, str) and is_simid(k))]
    if invalid_keys:
        msg = (
            f"invalid simid(s) found: {', '.join(repr(k) for k in invalid_keys)}. "
            r"simids must match the pattern [-\w]+ "
            "(letters, digits, underscores, hyphens only — dots are forbidden)"
        )
        raise SimflowConfigError(msg, block)


def query_runlist_db(metadata: MetadataRepository, query: str) -> list[str]:
    """Query the runlist DB stored in legend-datasets.

    Run expressions of the form ``r00n..r00m`` are automatically expanded into
    full run lists. If for example ``metadata.datasets.runlists.valid.phy.p02
    == "r000..r002"``:

    >>> query_runlist_db(metadata, "valid.phy.p02")
    ["l200-p02-r000-phy", "l200-p02-r001-phy", "l200-p02-r002-phy"]

    Parameters
    ----------
    metadata
        LEGEND metadata instance.
    query
        expression in the form `<tag>.<datatype>.<period>` (see contents of
        ``runlists.yaml`` in legend-datasets.

    """
    group, dtype, period = re.split(r"\W+", query)

    run_exprs = metadata.datasets.runlists[group][dtype][period]

    if not isinstance(run_exprs, list | tuple):
        run_exprs = [run_exprs]

    runs: list[str] = []
    for item in run_exprs:
        m = re.match(r"^r(\d+)\.\.r(\d+)$", item)
        if m is not None:
            r1, r2 = m.groups()
            runs.extend(
                [f"l200-{period}-r{r:03d}-{dtype}" for r in range(int(r1), int(r2) + 1)]
            )
        else:
            runs.append(item)

    return sorted(runs)


def expand_runlist(
    metadata: MetadataRepository, runlist: str | Iterable[str]
) -> list[str]:
    """Expands a runlist as passed to the Simflow configuration.

    A runlist is a list of:

    - runids in the form accepted by :func:`is_runid`;
    - runlist DB queries in the form ``<tag>.<datatype>.<period>`` (see
      :func:`query_runlist_db`).

    The returned list is sorted and free of duplicates.
    """
    if not isinstance(runlist, list | tuple):
        runlist = [runlist]

    runs = []
    for item in runlist:
        if item.startswith("~runlists:"):
            runs.extend(query_runlist_db(metadata, item.partition("~runlists:")[2]))
        else:
            if not is_runid(item):
                msg = f"{item} is not a valid runid"
                raise ValueError(msg)

            runs.append(item)

    # deduplicate: the same runid can be listed explicitly and also be part of a
    # runlist DB query. duplicates would silently corrupt the event partitioning
    # downstream (see partitioning.partition_simstat)
    return sorted(set(runs))


def expand_simlist(config: SimflowConfig, simlist: str | Iterable[str]) -> list[str]:
    """Expand a simlist as passed to the Simflow configuration.

    A simlist is a list (or a comma-separated string) of:

    - items in the form ``<tier>.<simid>``;
    - ``~simlists:NAME``, replaced by the list stored under ``NAME`` for the
      experiment in ``simprod/config/simlists.yaml`` of the Simflow metadata.

    The order is preserved and duplicates are removed.
    """
    if isinstance(simlist, str):
        simlist = simlist.split(",")

    items: list[str] = []
    for raw in simlist:
        item = raw.strip()
        if item.startswith("~simlists:"):
            name = item.partition("~simlists:")[2].strip()
            try:
                named = config.metadata.simprod.config.simlists[config.experiment][name]
            except KeyError as e:
                msg = f"no simlist named {name!r} for experiment {config.experiment!r}"
                raise SimflowConfigError(msg, "simflow-config.simlist") from e
            items.extend(i.strip() for i in named)
        elif item:
            items.append(item)

    return list(dict.fromkeys(items))


def get_runlist(config: SimflowConfig, simid: str) -> list[str]:
    """Gets the runlist assigned to a simulation.

    If not overridden in the hit-tier `simconfig`, returns the global runlist
    stored in ``config.runlist``.
    """
    default = config.runlist
    simcfg = get_simconfig(config, "hit")

    try:
        simcfg = simcfg[simid]
        runlist = simcfg["runlist"]
    except KeyError as e:
        if default is not None:
            runlist = default
        else:
            key = e.args[0]
            path = f"simprod.config.tier.hit.{config.experiment}.simconfig"
            if key == "runlist":
                path += f".{simid}"

            msg = f"'{key}' key not found and config.runlist fallback undefined"
            raise SimflowConfigError(msg, path) from e

    return expand_runlist(config.metadata, runlist)


# FIXME: this should be removed once the PRL25 data is reprocessed
def _get_lh5_table(
    metadata: MetadataRepository,
    fname: str | Path,
    hpge: str,
    tier: str,
    runid: str,
    use_rawid: bool = False,
) -> str:
    """The correct LH5 table path.

    Determines the correct path to a `hpge` detector table in tier `tier`. With
    `use_rawid`, always return the old ``ch{rawid}/{tier}`` path.
    """
    # check if the latest format is available
    path = f"{tier}/{hpge}"
    if not use_rawid and lh5.ls(fname, path) == [path]:
        return path

    # otherwise fall back to the old format
    timestamp = runinfo(metadata, runid).start_key

    chmap = metadata.channelmap(timestamp)

    rawid = chmap[hpge].daq.rawid
    return f"ch{rawid}/{tier}"
