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
import json

from dbetto import AttrsDict
from snakemake.io import Wildcards

from . import SimflowConfig
from .exceptions import SimflowConfigError


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


def hash_dict(d):
    """Compute the hash of a Python dict."""
    if isinstance(d, AttrsDict):
        d = d.to_dict()

    s = json.dumps(d, sort_keys=True)
    return hashlib.sha256(s.encode()).hexdigest()


def smk_hash_simconfig(
    config: SimflowConfig,
    wildcards: Wildcards,
    field: str | None = None,
    ignore: list | None = None,
    **kwargs,
):
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

    return hash_dict(scfg)
