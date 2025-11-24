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
from __future__ import annotations

import shlex
from pathlib import Path

import legenddataflowscripts as lds
import numpy as np

from . import SimflowConfig, patterns
from .exceptions import SimflowConfigError
from .metadata import get_simconfig


def remage_run(
    config: SimflowConfig,
    simid: str,
    tier: str = "stp",
    geom: str | Path = "{input.geom}",
    procs: int = 1,
    output: str | Path = "{output}",
    macro_free: bool = False,
) -> str:
    """Build a remage CLI invocation string for a given simulation.

    This constructs a shell-escaped command line for remage by first rendering
    the macro via :func:`make_remage_macro` using the simulation configuration
    (from ``simconfig.yaml``), and then assembling the remage CLI with the
    appropriate arguments and macro handling.

    Notes
    -----
    - Compatible with remage >= v0.16.
    - When ``macro_free`` is False (default), the command passes the macro file
      path and supplies macro substitutions via ``--macro-substitutions``.
    - When ``macro_free`` is True, the rendered macro content is inlined on the
      CLI (comments and empty lines removed) and values are pre-substituted.
    - Two substitutions are always provided:
      ``N_EVENTS`` (from ``primaries_per_job`` or benchmark override) and
      ``SEED`` (a random 32-bit integer).
    - If ``config.runcmd.remage`` is set, it is used to determine the remage
      executable (split with :func:`shlex.split`), otherwise ``remage`` is used.

    Parameters
    ----------
    config
        Snakemake-like configuration mapping. Must include metadata required by
        :func:`make_remage_macro` and optional ``benchmark`` and ``runcmd``
        sections.
    simid
        Simulation identifier for which to construct the command.
    tier
        Simulation tier (e.g., ``"stp"``, ``"ver"``). Default is ``"stp"``.
    geom
        Path (or Snakemake placeholder) to the GDML geometry file.
    procs
        Number of threads to pass to remage (integer or Snakemake placeholder).
        Internally uses remage's ``--procs``.
    output
        Path (or Snakemake placeholder) to the output remage file.
    macro_free
        If True, inline the macro contents on the CLI; if False, reference the
        macro file and pass substitutions via ``--macro-substitutions``.

    Returns
    -------
    A shell-escaped command line suitable for direct execution.
    """

    # get the config block for this tier/simid
    block = f"simprod.config.tier.{tier}.{config.experiment}.simconfig.{simid}"
    sim_cfg = get_simconfig(config, tier, simid=simid)

    # get macro
    macro_text, _ = make_remage_macro(config, simid, tier=tier)

    # need some modifications if this is a benchmark run
    try:
        n_prim_pj = sim_cfg.primaries_per_job
        is_benchmark = False

        if "benchmark" in config:
            is_benchmark = config.benchmark.get("enabled", False)
            if is_benchmark:
                n_prim_pj = config.benchmark.n_primaries[tier]
    except KeyError as e:
        msg = f"key {e} not found!"
        raise SimflowConfigError(msg, block) from e

    # substitution rules
    cli_subs = {
        "N_EVENTS": int(n_prim_pj / int(procs)),
        # TODO: check correct range
        "SEED": np.random.default_rng().integers(
            0, np.iinfo(np.int32).max + 1, dtype=np.uint32
        ),
    }

    # prepare command line
    if "runcmd" in config and "remage" in config.runcmd:
        remage_exe = shlex.split(config.runcmd.get("remage", "remage"))
    else:
        remage_exe = ["remage"]

    cmd = [
        *remage_exe,
        "--ignore-warnings",
        "--merge-output-files",
        "--log-level=detail",
        "--procs",
        str(procs),
        "--gdml-files",
        str(geom),
        "--output-file",
        str(output),
    ]

    if macro_free:
        # actually substitute values in the command line
        for k, v in cli_subs.items():
            macro_text = macro_text.replace(f"{{{k}}}", str(v))

        cmd += ["--"]

        # remove empty lines and comments if passing macro commands to the cli
        cmd += [
            line for line in macro_text.splitlines() if line.strip() and line[0] != "#"
        ]
    else:
        # otherwise just pass the macro file
        cmd += ["--macro-substitutions"]
        cmd += [shlex.quote(f"{k}={v}") for k, v in cli_subs.items()]

        cmd += ["--"]

        cmd += [
            patterns.input_simjob_filename(config, tier=tier, simid=simid).as_posix()
        ]

    return shlex.join(cmd)


def make_remage_macro(
    config: SimflowConfig, simid: str, tier: str = "stp"
) -> (str, Path):
    """Render the remage macro for a given simulation and write it to disk.

    This function reads the simulation configuration for the provided tier/simid,
    assembles the macro substitutions (e.g. ``GENERATOR``, ``CONFINEMENT``)
    using values and references defined under config.metadata, renders the
    specified macro template, writes the final macro file to the canonical
    input path, and returns both the macro text and the output file path.

    Parameters
    ----------
    config
        Mapping-like Snakemake configuration that supports attribute-style access
        (e.g. ``config.experiment``, ``config.metadata``, etc.). The following fields
        are used:
        - ``experiment``: name of the experiment to select tier-specific metadata.
        - ``metadata.tier[tier][experiment].generators``: generator definitions.
        - ``metadata.tier[tier][experiment].confinement``: confinement definitions.
    simid
        Simulation identifier to select the simconfig.
    tier
        Simulation tier (e.g. "stp", "ver", ...). Default is "stp".

    Returns
    -------
    A tuple with:
    - The rendered macro text.
    - The path where the macro was written.

    Notes
    -----
    - The macro template path is taken from the simconfig `template` field.
    - Supported substitutions currently include: ``GENERATOR`` and
      ``CONFINEMENT``.
    - The user can provide arbitrary macro substitutions with the
      optional `macro_substitutions` field.
    - The macro is written to the canonical path returned by
      :func:`.patterns.input_simjob_filename`.
    """
    # get the config block for this tier/simid
    block = f"simprod.config.tier.{tier}.{config.experiment}.simconfig.{simid}"
    sim_cfg = get_simconfig(config, tier, simid=simid)
    mac_subs = {}

    # determine whether external vertices are required
    if "vertices" in sim_cfg:
        raise NotImplementedError()
        mac_subs.update({"VERTICES_FILE": None})

    # configure generator
    if "generator" in sim_cfg:
        if not isinstance(sim_cfg.generator, str) or not sim_cfg.generator.startswith(
            "~defines:"
        ):
            msg = (
                "the field must be a string prefixed by ~define:",
                f"{block}.generator",
            )
            raise SimflowConfigError(*msg)

        key = sim_cfg.generator.removeprefix("~defines:")
        try:
            generator = config.metadata.simprod.config.tier[tier][
                config.experiment
            ].generators[key]
        except KeyError as e:
            msg = f"key {e} not found!"
            raise SimflowConfigError(msg, block) from e

        if not isinstance(generator, str):
            generator = "\n".join(generator)

        mac_subs["GENERATOR"] = generator

    # configure generator
    if "confinement" in sim_cfg:
        confinement = None
        if isinstance(sim_cfg.confinement, str):
            if sim_cfg.confinement.startswith("~defines:"):
                key = sim_cfg.confinement.removeprefix("~defines:")
                try:
                    confinement = config.metadata.simprod.config.tier[tier][
                        config.experiment
                    ].confinement[key]
                except KeyError as e:
                    msg = f"key {e} not found!"
                    raise SimflowConfigError(msg, block) from e

            elif sim_cfg.confinement.startswith(
                ("~volumes.surface:", "~volumes.bulk:")
            ):
                confinement = [
                    "/RMG/Generator/Confine Volume",
                    "/RMG/Generator/Confinement/Physical/AddVolume "
                    + sim_cfg.confinement.partition(":")[2],
                ]
                if sim_cfg.confinement.startswith("~volumes.surface:"):
                    confinement += ["/RMG/Generator/Confinement/SampleOnSurface true"]
            else:
                confinement = None

        elif isinstance(sim_cfg.confinement, list | tuple):
            confinement = ["/RMG/Generator/Confine Volume"]
            for val in sim_cfg.confinement:
                if val.startswith(("~volumes.surface:", "~volumes.bulk:")):
                    confinement += [
                        "/RMG/Generator/Confinement/Physical/AddVolume "
                        + val.partition(":")[2]
                    ]
                    if val.startswith("~volumes.surface:"):
                        confinement += [
                            "/RMG/Generator/Confinement/SampleOnSurface true"
                        ]
                else:
                    confinement = None

        if confinement is None:
            msg = (
                (
                    "the field must be a str or list[str] prefixed by "
                    "~define: / ~volumes.surface: / ~volumes.bulk:"
                ),
                f"{block}.confinement",
            )
            raise SimflowConfigError(*msg)

        if not isinstance(confinement, str):
            confinement = "\n".join(confinement)

        mac_subs["CONFINEMENT"] = confinement

    # the user might want to substitute some other variables
    mac_subs |= sim_cfg.get("macro_substitutions", {})

    # read in template and substitute
    template_path = get_simconfig(config, tier, simid=simid, field="template")
    with Path(template_path).open() as f:
        text = lds.subst_vars(f.read().strip(), mac_subs, ignore_missing=False)

    # now write the macro to disk
    ofile = patterns.input_simjob_filename(config, tier=tier, simid=simid)
    ofile.parent.mkdir(parents=True, exist_ok=True)
    with ofile.open("w") as f:
        f.write(text)

    return text, ofile
