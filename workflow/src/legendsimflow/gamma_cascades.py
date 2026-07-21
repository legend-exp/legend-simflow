# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>
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
"""Build the _remage_ macro block registering the MAURINA gamma cascades.

The ``maurina_gamma_cascades`` field of the `stp` tier simconfig enables the
tabulated MAURINA gamma cascades: they change how the compound nucleus
de-excites after a neutron capture. It can be either ``true`` (register all
isotopes) or a list of ``(Z, A)`` pairs selecting a subset.

The commands must be issued **before** ``/run/initialize``, since they are read
when ``RMGPhysics::ConstructProcess()`` builds the high-precision hadronic
models, which Geant4 does at initialization time.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from . import SimflowConfig, nersc
from .exceptions import SimflowConfigError

#: Glob matching the per-isotope filelists in the cascade data repository. The
#: basenames are not uniform across isotopes (``ge76_``, ``gd155_``, ``cr53_``,
#: ...), so they cannot be constructed from `Z` and `A`.
FILELIST_GLOB = "*/*/*_ncapture_filelist.txt"

#: Name of the macro template variable holding the rendered block.
MACRO_VARIABLE = "MAURINA_GAMMA_CASCADES"


def maurina_macro_commands(
    config: SimflowConfig,
    block: str | None = None,
    subset: set[tuple[int, int]] | None = None,
) -> list[str]:
    """Return the commands registering all MAURINA gamma cascades.

    Walks the ``<Z>/<A>/`` directory structure of the cascade data repository at
    ``config.paths.maurina_gamma_cascades`` and registers the filelist found in
    each.

    Note
    ----
    The filelists are referenced at their location in the repository, because
    _remage_ resolves the cascade file names they contain relative to the
    directory holding the filelist.
    """
    if "maurina_gamma_cascades" not in config.paths:
        msg = (
            "MAURINA gamma cascades are enabled but "
            "'paths.maurina_gamma_cascades' is not configured: it must point to "
            "a checkout of the legend-exp/generated-maurina-output repository"
        )
        raise SimflowConfigError(msg, block)

    root = Path(nersc.dvs_ro(config, config.paths.maurina_gamma_cascades))
    if not root.is_dir():
        msg = f"MAURINA gamma cascade data directory '{root}' does not exist"
        raise SimflowConfigError(msg, block)

    filelists = sorted(root.glob(FILELIST_GLOB))
    if not filelists:
        msg = (
            f"no '{FILELIST_GLOB}' found below MAURINA gamma cascade directory '{root}'"
        )
        raise SimflowConfigError(msg, block)

    commands = ["/RMG/Processes/UseGrabmayrsGammaCascades true"]
    seen: dict[tuple[int, int], Path] = {}

    for filelist in filelists:
        z_dir, a_dir = filelist.parent.parent.name, filelist.parent.name
        try:
            isotope = (int(z_dir), int(a_dir))
        except ValueError as e:
            msg = (
                f"cannot parse Z/A from MAURINA gamma cascade path '{filelist}': "
                f"the two directory levels below '{root}' must be integers, got "
                f"'{z_dir}/{a_dir}'"
            )
            raise SimflowConfigError(msg, block) from e

        # ignore isotopes when a subset is defined and the isotope is not contained in it
        if subset is not None and isotope not in subset:
            continue

        # remage resets previously registered cascades for an isotope, so a
        # second filelist would silently shadow the first
        if isotope in seen:
            msg = (
                f"multiple MAURINA gamma cascade filelists for isotope "
                f"Z={isotope[0]} "
                f"A={isotope[1]}: '{seen[isotope]}' and '{filelist}'"
            )
            raise SimflowConfigError(msg, block)
        seen[isotope] = filelist

        commands.append(
            "/RMG/GrabmayrGammaCascades/SetGammaCascadeFilelist "
            f"{isotope[0]} {isotope[1]} {filelist}"
        )

    commands.append("/RMG/GrabmayrGammaCascades/SetGammaCascadeRandomStartLocation 1")
    return commands


def maurina_macro_block(
    config: SimflowConfig, sim_cfg: Mapping, block: str | None = None
) -> str:
    """Return the full macro block for a simulation, or ``""`` if unconfigured.

    Renders :func:`maurina_macro_commands` if the optional
    ``maurina_gamma_cascades`` field of `sim_cfg` is set.
    """
    if "maurina_gamma_cascades" not in sim_cfg:
        return ""

    if isinstance(  # test for bool
        sim_cfg["maurina_gamma_cascades"], bool
    ):
        if not sim_cfg.get("maurina_gamma_cascades", False):
            return ""
        return "\n".join(maurina_macro_commands(config, block))

    if isinstance(  # test for list
        sim_cfg["maurina_gamma_cascades"], list
    ):
        subset = set()
        for pair in sim_cfg["maurina_gamma_cascades"]:
            if (
                not isinstance(pair, (list, tuple))
                or len(pair) != 2
                or not all(isinstance(x, int) for x in pair)
            ):
                msg = (
                    "the 'maurina_gamma_cascades' list must contain pairs of integers (Z, A), "
                    f"got {pair!r}"
                )
                raise SimflowConfigError(msg, block)

            z, a = pair
            subset.add((z, a))
        return "\n".join(maurina_macro_commands(config, block, subset=subset))

    # neither bool nor list
    msg = (
        "the 'maurina_gamma_cascades' field must be either a boolean or a "
        f"list of (Z, A) pairs, got {sim_cfg['maurina_gamma_cascades']}"
    )
    raise SimflowConfigError(msg, block)
