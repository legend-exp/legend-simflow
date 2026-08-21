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
from legendmeta import LegendMetadata
from legendmeta.police import validate_dict_schema
from snakemake.iocontainers import Wildcards

from . import SimflowConfig, utils
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

#: runid prefix assumed for an experiment whose identifier does not carry one.
#: See :func:`experiment_prefix`.
DEFAULT_RUNID_PREFIX = "l200"

#: first letter of a detector name, by detector type. See
#: :func:`get_crystal_name`.
CRYSTAL_TYPE_IDS = {"bege": "B", "coax": "C", "ppc": "P", "icpc": "V"}

#: subdirectory of ``paths.config`` that holds the metadata overlay of each
#: experiment. See :func:`load_metadata_overlay`.
METADATA_OVERLAY_DIR = "metadata"

#: errors a metadata database raises when it does not hold the queried item.
#: :meth:`dbetto.catalog.Catalog.valid_for` raises ``RuntimeError`` for a
#: timestamp that precedes every validity entry of a directory. This is how a
#: channel map lookup misses. See :func:`lookup`.
LOOKUP_ERRORS = (LookupError, FileNotFoundError, RuntimeError)


def metadata_overlay_dirname(config: SimflowConfig) -> Path:
    """The metadata overlay directory of the current experiment."""
    return Path(config.paths.config) / METADATA_OVERLAY_DIR / config.experiment


def load_metadata_overlay(
    config: SimflowConfig, *, logger: logging.Logger | None = None
) -> LegendMetadata | None:
    """Load the metadata overlay of the current experiment.

    An experiment that does not exist yet is absent from `legend-metadata`. Its
    metadata lives in `legend-simflow-config` instead, under
    ``{paths.config}/metadata/{experiment}/``. The tree uses the
    `legend-metadata` layout, so the same class reads it.

    The Simflow queries `legend-metadata` first and the overlay second (see
    :func:`lookup`). An experiment without such a directory never reaches the
    overlay.

    Returns ``None`` if the directory does not exist.

    Parameters
    ----------
    config
        Simflow configuration object, with ``config.metadata`` already attached.
    logger
        Logger to use for status messages (e.g. the Snakemake logger when called
        from a Snakefile). Defaults to the module logger.

    """
    log_ = logger if logger is not None else log

    if "config" not in config.paths:
        return None

    path = metadata_overlay_dirname(config)
    if not path.is_dir():
        return None

    msg = f"loading the metadata overlay of {config.experiment} from {path}"
    log_.info(msg)

    overlay = LegendMetadata(path)
    validate_metadata_overlay(config.metadata, overlay, path)
    return overlay


def validate_metadata_overlay(
    metadata: LegendMetadata, overlay: LegendMetadata, path: Path
) -> None:
    """Refuse an overlay whose runs collide with the `legend-metadata` runs.

    The Simflow looks a run up by period and run number alone. A period that
    `legend-metadata` also defines therefore resolves there first, and gives the
    wrong run. :func:`lookup` cannot detect this collision. An overlay must
    number its periods outside the range of the real experiment.

    Parameters
    ----------
    metadata
        LEGEND metadata database.
    overlay
        Metadata overlay of the experiment.
    path
        Directory that holds `overlay`. The error message names it.

    """
    try:
        periods = set(overlay.datasets.runinfo)
    except LOOKUP_ERRORS:
        return

    try:
        existing = set(metadata.datasets.runinfo)
    except LOOKUP_ERRORS:
        return

    clash = sorted(periods & existing)
    if clash:
        msg = (
            f"legend-metadata also defines the period(s) {', '.join(clash)} of "
            f"the metadata overlay in {path}. A run of such a period resolves "
            "to the legend-metadata run. Number the overlay periods outside the "
            "range of the real experiment, for example p99"
        )
        raise SimflowConfigError(msg, "paths.config")


def lookup(
    config: SimflowConfig,
    query: Callable[[LegendMetadata], object],
    default: object = _MISSING,
) -> object:
    """Run `query` against the metadata databases, in order.

    The function queries `legend-metadata` first and the metadata overlay second
    (see :func:`load_metadata_overlay`). The real database therefore always
    wins.

    A database that raises one of :data:`LOOKUP_ERRORS` does not hold the item,
    so the function tries the next one. If no database holds the item, the
    function raises the error of the last one.

    Parameters
    ----------
    config
        Simflow configuration object.
    query
        Callable that takes a metadata database and returns the queried item.
    default
        Value to return when no database holds the item. Without it, the
        function raises.

    """
    dbs = [config.metadata]

    overlay = config.get("metadata_overlay")
    if overlay is not None:
        dbs.append(overlay)

    error = None
    for db in dbs:
        try:
            return query(db)
        except LOOKUP_ERRORS as e:
            error = e

    if default is not _MISSING:
        return default

    raise error


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
        # an empty simconfig.yaml is a valid configuration (e.g. a hit tier that
        # takes its runlist from config.runlist), but YAML loads it as None
        simcfg = _m.tier[tier][config.experiment].simconfig or AttrsDict({})

        if simid is None:
            block = f"simprod.config.tier.{tier}.{config.experiment}"
            validate_simconfig_keys(simcfg, block + ".simconfig")
            return simcfg
        if field is None:
            return simcfg[simid]
        return simcfg[simid][field]

    except KeyError as e:
        msg = f"key {e} not found!"
        raise SimflowConfigError(msg, block) from e
    except FileNotFoundError as e:
        raise SimflowConfigError(e, block) from e


def get_tier_settings(config: SimflowConfig, tier: str) -> AttrsDict:
    """Return the settings block for *tier* and the current experiment."""
    return config.metadata.simprod.config.tier[tier][config.experiment].settings


def get_par_settings(config: SimflowConfig, par: str) -> AttrsDict:
    """Return the settings block for *par* and the current experiment."""
    if (
        config.experiment in config.metadata.simprod.config.pars
        and par in config.metadata.simprod.config.pars[config.experiment].geds
        and "settings"
        in config.metadata.simprod.config.pars[config.experiment].geds[par]
    ):
        return config.metadata.simprod.config.pars[config.experiment].geds[par].settings
    return AttrsDict({})


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
    config: SimflowConfig, det_name: str, runid: str, default: str | None = None
) -> str:
    """Get the usability for analysis of `det_name` in run `runid`.

    Looks for the ``analysis.usability`` metadata field in the channel map. By
    default, an error is thrown if no information is found. If `default` is set
    to a non-None value, it will be returned.
    """
    rinfo = get_runinfo(config, runid)
    chmap = get_channelmap(config, rinfo.start_key)
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


def get_runinfo(config: SimflowConfig, runid: str) -> AttrsDict:
    """Get the `datasets.runinfo` entry for a LEGEND run identifier.

    Parameters
    ----------
    config
        Simflow configuration object.
    runid
        a run identifier in the format ``<experiment>-<period>-<run>-<datatype>``.

    """
    _, period, run, datatype = re.split(r"\W+", runid)
    return lookup(config, lambda db: db.datasets.runinfo[period][run][datatype])


def get_channelmap(config: SimflowConfig, timestamp: str) -> AttrsDict:
    """The channel map of the current experiment, valid at `timestamp`.

    Queries the metadata overlay when `legend-metadata` does not hold the
    entry. See :func:`lookup`.
    """
    return lookup(config, lambda db: db.channelmap(timestamp, skip_version_check=True))


def get_diode(
    config: SimflowConfig, det_name: str, default: object = _MISSING
) -> AttrsDict:
    """The `hardware.detectors.germanium.diodes` entry of `det_name`.

    Queries the metadata overlay when `legend-metadata` does not hold the
    entry. See :func:`lookup`.
    """
    return lookup(
        config,
        lambda db: db.hardware.detectors.germanium.diodes[det_name],
        default=default,
    )


def get_crystal(
    config: SimflowConfig, crystal_name: str, default: object = _MISSING
) -> AttrsDict:
    """The `hardware.detectors.germanium.crystals` entry of `crystal_name`.

    Queries the metadata overlay when `legend-metadata` does not hold the
    entry. See :func:`lookup`.
    """
    return lookup(
        config,
        lambda db: db.hardware.detectors.germanium.crystals[crystal_name],
        default=default,
    )


def reference_cal_run(config: SimflowConfig, runid: str) -> str:
    """The reference calibration run for `runid`.

    Warning
    -------
    This function does not account for dataflow overrides (e.g. calibration
    back-applying)!

    """
    exp, period, run, datatype = re.split(r"\W+", runid)

    msg = f"looking for reference calibration run for {runid}"
    log.info(msg)

    # TODO: this is a temporary workaround for the fact that p16-r007 and later runs should use p16-r006
    if period == "p16" and (run > "r006"):
        run = "r006"
        runid = f"{exp}-{period}-{run}-{datatype}"
        log.info("p16-r006 + so we set the cal run to r006")

    if datatype == "cal":
        return runid

    p_runinfo = lookup(config, lambda db: db.datasets.runinfo[period])

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
    config: SimflowConfig,
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
    config
        Simflow configuration object.
    par
        name of directory under ``metadata.simprod.config.pars.{experiment}``.
        Can be a nested property, as in e.g. ``geds.opv.value``. ``.`` and
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
        directory = config.metadata["simprod/config/pars"][experiment][par]
        return directory.on(get_runinfo(config, runid).start_key, category=datatype)
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


def get_sanitized_fccd(config: SimflowConfig, det_name: str) -> float:
    """Return the FCCD value for `det_name`, falling back to 1 mm if the FCCD field is absent.

    Parameters
    ----------
    config
        Simflow configuration object.
    det_name
        Detector name.

    """
    det_meta = get_diode(config, det_name, default=None)
    if det_meta is None:
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


def experiment_prefix(experiment: str) -> str:
    """Get the runid prefix of an experiment.

    The prefix is the leading letters-then-digits part of the experiment
    identifier, i.e. the name of the setup itself, stripped of the
    configuration tag: ``l200cfg09`` and ``l1000dsg01`` give ``l200`` and
    ``l1000``. Falls back to :data:`DEFAULT_RUNID_PREFIX` for identifiers that
    do not follow the convention.

    Parameters
    ----------
    experiment
        Experiment identifier (e.g. ``l200cfg09``, ``l1000dsg01``).

    """
    match = re.match(r"^[A-Za-z]+\d+", experiment)
    if match is None:
        msg = f"cannot determine the runid prefix of experiment {experiment}, assuming {DEFAULT_RUNID_PREFIX}"
        log.warning(msg)
        return DEFAULT_RUNID_PREFIX

    return match.group(0)


def query_runlist_db(
    metadata: LegendMetadata, query: str, prefix: str = DEFAULT_RUNID_PREFIX
) -> list[str]:
    """Query the runlist DB stored in legend-datasets.

    Run expressions of the form ``r00n..r00m`` are automatically expanded into
    full run lists. If for example ``metadata.datasets.runlists.valid.phy.p02
    == "r000..r002"``:

    >>> query_runlist_db(metadata, "valid.phy.p02", "l200")
    ["l200-p02-r000-phy", "l200-p02-r001-phy", "l200-p02-r002-phy"]

    Parameters
    ----------
    metadata
        LEGEND metadata instance.
    query
        expression in the form `<tag>.<datatype>.<period>` (see contents of
        ``runlists.yaml`` in legend-datasets.
    prefix
        Runid prefix, see :func:`experiment_prefix`.

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
                [
                    f"{prefix}-{period}-r{r:03d}-{dtype}"
                    for r in range(int(r1), int(r2) + 1)
                ]
            )
        else:
            runs.append(item)

    return sorted(runs)


def expand_runlist(
    metadata: LegendMetadata,
    runlist: str | Iterable[str],
    prefix: str = DEFAULT_RUNID_PREFIX,
) -> list[str]:
    """Expands a runlist as passed to the Simflow configuration.

    A runlist is a list of:

    - runids in the form accepted by :func:`is_runid`;
    - runlist DB queries in the form ``<tag>.<datatype>.<period>`` (see
      :func:`query_runlist_db`).

    Parameters
    ----------
    metadata
        LEGEND metadata instance.
    runlist
        The runlist to expand.
    prefix
        Runid prefix, see :func:`experiment_prefix`.

    """
    if not isinstance(runlist, list | tuple):
        runlist = [runlist]

    runs = []
    for item in runlist:
        if item.startswith("~runlists:"):
            runs.extend(
                query_runlist_db(metadata, item.partition("~runlists:")[2], prefix)
            )
        else:
            if not is_runid(item):
                msg = f"{item} is not a valid runid"
                raise ValueError(msg)

            runs.append(item)
    return sorted(runs)


def get_runlist(config: SimflowConfig, simid: str) -> list[str]:
    """Gets the runlist assigned to a simulation.

    If not overridden in the hit-tier `simconfig`, returns the global runlist
    stored in ``config.runlist``.
    """
    default = config.get("runlist", None)
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

    prefix = experiment_prefix(config.experiment)
    return lookup(config, lambda db: expand_runlist(db, runlist, prefix=prefix))


# FIXME: this should be removed once the PRL25 data is reprocessed
def _get_lh5_table(
    config: SimflowConfig,
    fname: str | Path,
    hpge: str,
    tier: str,
    runid: str,
) -> str:
    """The correct LH5 table path.

    Determines the correct path to a `hpge` detector table in tier `tier`.
    """
    # check if the latest format is available
    path = f"{tier}/{hpge}"
    if lh5.ls(fname, path) == [path]:
        return path

    # otherwise fall back to the old format
    chmap = get_channelmap(config, get_runinfo(config, runid).start_key)

    rawid = chmap[hpge].daq.rawid
    return f"ch{rawid}/{tier}"


def get_crystal_name(diode_meta: AttrsDict) -> str:
    """Get the name of the crystal an HPGe detector was cut from.

    Assembled from the detector type and the crystal production information: the
    detector ``V05261B`` (an ICPC, order 5, crystal 261) was cut from crystal
    ``V05261``.

    Parameters
    ----------
    diode_meta
        Diode metadata, i.e. an entry of
        ``hardware.detectors.germanium.diodes``.

    """
    return (
        CRYSTAL_TYPE_IDS[diode_meta.type]
        + format(diode_meta.production.order, "02d")
        + diode_meta.production.crystal
    )
