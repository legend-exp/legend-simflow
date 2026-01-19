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
from collections.abc import Collection
from pathlib import Path

from dbetto import AttrsDict
from legendmeta import LegendMetadata
from legendmeta.police import validate_dict_schema
from lgdo import lh5
from snakemake.io import Wildcards

from . import SimflowConfig, utils
from .exceptions import SimflowConfigError

log = logging.getLogger(__name__)


def get_simconfig(
    config: SimflowConfig,
    tier: str,
    simid: str | None = None,
    field: str | None = None,
) -> AttrsDict:
    """Get the simulation configuration.

    Gets the simconfig and throws proper exceptions.

    Parameters
    ----------
    config
        Snakemake config.
    tier
        tier name.
    simid
        simulation identifier.
    field
        if not none, return the value of this key in the simconfig.
    """
    try:
        _m = config.metadata.simprod.config
    except FileNotFoundError as e:
        raise SimflowConfigError(e) from e

    block = f"simprod.config.tier.{tier}.{config.experiment}.simconfig"
    try:
        if simid is None:
            block = f"simprod.config.tier.{tier}.{config.experiment}"
            return _m.tier[tier][config.experiment].simconfig
        if field is None:
            return _m.tier[tier][config.experiment].simconfig[simid]
        return _m.tier[tier][config.experiment].simconfig[simid][field]

    except KeyError as e:
        msg = f"key {e} not found!"
        raise SimflowConfigError(msg, block) from e
    except FileNotFoundError as e:
        raise SimflowConfigError(e, block) from e


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
        if not none, return the value of this key in the simconfig.
    ignore
        exclude these fields from the hash.
    kwargs
        provide a value for wildcards that might not be present in `wildcards`.
    """
    tier = kwargs["tier"] if "tier" in kwargs else wildcards.tier  # noqa: SIM401
    simid = kwargs["simid"] if "simid" in kwargs else wildcards.simid  # noqa: SIM401

    scfg = get_simconfig(config, tier, simid)

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
    with file_path.open() as f:
        return int(f.read().strip())


def usability(
    metadata: LegendMetadata, det_name: str, runid: str, default: str | None = None
) -> str:
    """Get the usability for analysis of `det_name` in run `ruinid`.

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


def runinfo(metadata: LegendMetadata, runid: str) -> str:
    """Get the `datasets.runinfo` entry for a LEGEND run identifier.

    Parameters
    ----------
    metadata
        LEGEND metadata database.
    runid
        a run identifier in the format ``<experiment>-<period>-<run>-<datatype>``.
    """
    _, period, run, datatype = re.split(r"\W+", runid)
    return metadata.datasets.runinfo[period][run][datatype]


def reference_cal_run(metadata: LegendMetadata, runid: str) -> str:
    """The reference calibration run for `runid`.

    Warning
    -------
    This function does not account for dataflow overrides (e.g. calibration
    back-applying)!
    """
    exp, period, run, datatype = re.split(r"\W+", runid)

    if datatype == "cal":
        return runid

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


def simpars(metadata: LegendMetadata, par: str, runid: str) -> AttrsDict:
    """Extract simflow parameters for a certain LEGEND run.

    Queries the simflow parameters database stored under
    ``simprod.config.pars`` by parameter name `par` and LEGEND run identifier
    `runid`.

    Parameters
    ----------
    metadata
        LEGEND metadata database.
    par
        name of directory under ``metadata.simprod.config.pars``. Can be a
        nested property, as in e.g. ``geds.opv.value``. ``.`` and ``/`` are
        allowed separators.
    runid
        a run identifier in the format ``<experiment>-<period>-<run>-<datatype>``.
    """
    par = par.replace(".", "/")
    datatype = re.split(r"\W+", runid)[-1]
    directory = metadata["simprod/config/pars"][par]
    return directory.on(runinfo(metadata, runid).start_key, system=datatype)


def get_vtx_simconfig(config, simid):
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


def get_sanitized_fccd(metadata: LegendMetadata, det_name: str) -> float:
    det_meta = metadata.hardware.detectors.germanium.diodes[det_name]

    has_fccd_meta = validate_dict_schema(
        det_meta.characterization,
        {"combined_0vbb_analysis": {"fccd_in_mm": 0}},
        greedy=False,
        verbose=False,
    )

    if not has_fccd_meta:
        msg = f"{det_name} metadata does not seem to contain usable FCCD data, setting to 1 mm"
        log.warning(msg)
        fccd = 1
    else:
        fccd = det_meta.characterization.combined_0vbb_analysis.fccd_in_mm

    return fccd


def is_runid(runid: str) -> bool:
    """Is a runid (run identifier) correctly formatted?

    It should be in the form
    ``l200-<period>-<run>-<datatype>``/``l200-pNN-rMMM-AAA``.
    """
    return re.match(r"^l200-p(\d{2})-r(\d{3})-([A-Za-z]+)$", runid) is not None


def query_runlist_db(metadata: LegendMetadata, query: str) -> list[str]:
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
    metadata: LegendMetadata, runlist: str | Collection[str]
) -> list[str]:
    """Expands a runlist as passed to the Simflow configuration.

    A runlist is a list of:

    - runids in the form accepted by :func:`is_runid`;
    - runlist DB queries in the form ``<tag>.<datatype>.<period>`` (see
      :func:`query_runlist_db`).
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

    return expand_runlist(config.metadata, runlist)


# FIXME: this should be removed once the PRL25 data is reprocessed
def _get_lh5_table(
    metadata: LegendMetadata,
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
    timestamp = runinfo(metadata, runid).start_key
    log.info(timestamp)

    chmap = metadata.channelmap(timestamp)
    log.info(chmap.keys())
    rawid = chmap[hpge].daq.rawid
    return f"ch{rawid}/{tier}"
