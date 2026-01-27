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
from copy import copy
from pathlib import Path

import legenddataflowscripts as lds

from . import SimflowConfig, nersc, patterns, utils
from .exceptions import SimflowConfigError
from .metadata import get_simconfig


def remage_run(
    config: SimflowConfig,
    simid: str,
    *,
    jobid: str | None = None,
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
    - ``SEED`` is meant to be used as remage seed. It is determined by
      converting `output` to a 32-bit integer hash. If provided, the user
      ``config.simflow_rng_seed`` integer is added as offset.
    - The ``JOBID`` substitution is also provided if the `jobid` argument is
      not ``None``.
    - If ``config.runcmd.remage`` is set, it is used to determine the remage
      executable (split with :func:`shlex.split`), otherwise ``remage`` is used.
    - If ``config.nersc.dvs_ro`` is set, remage is set to read all inputs from
      the read-only filesystem mount ``/dvs/ro`` at NERSC.
    - If ``config.nersc.scratch`` is set, the command will write the output
      file on the scratch disk and move it to the final expected destination at
      the end.

    Parameters
    ----------
    config
        Snakemake-like configuration mapping. Must include metadata required by
        :func:`make_remage_macro` and optional ``benchmark`` and ``runcmd``
        sections.
    simid
        Simulation identifier for which to construct the command.
    jobid
        Job identifier for the simulation run (string holding a zero-padded
        integer). Used as remage CLI macro substitution in case the macro
        contains it (e.g. if a vertices file is used).
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
    if not isinstance(output, Path):
        output = Path(str(output))

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

    # idea behind seed logic: we want the simflow to be reproducible across
    # runs. we convert the output file name to a unique integer and use that as
    # remage seed. we still want the user to still select some different global
    # seed, so we offset it with an integer optionally specified in the simflow
    # config

    # substitution rules
    cli_subs = {
        "N_EVENTS": int(n_prim_pj / int(procs)),
        "SEED": utils.string_to_remage_seed(
            str(output), seed=config.get("simflow_rng_seed", 0)
        ),
    }
    if jobid is not None:
        cli_subs["JOBID"] = jobid

    # prepare command line
    if "runcmd" in config and "remage" in config.runcmd:
        remage_exe = shlex.split(config.runcmd.get("remage", "remage"))
    else:
        remage_exe = ["remage"]

    output_final = copy(output)
    if nersc.is_scratch_enabled(config):
        output = nersc.on_scratch(config, output)

    cmd = [
        *remage_exe,
        "--ignore-warnings",
        "--merge-output-files",
        "--overwrite",
        "--log-level=detail",
        "--procs",
        str(procs),
        "--gdml-files",
        str(nersc.dvs_ro(config, geom)),
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

    joined_cmd = shlex.join(cmd)

    if nersc.is_scratch_enabled(config):
        joined_cmd += f" && mv -v {output} {output_final}"

    return joined_cmd


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
    - If ``config.nersc.dvs_ro`` is set, the vertices file will be read from
      the read-only filesystem mount ``/dvs_ro`` at NERSC.
    """
    # get the config block for this tier/simid
    block = f"simprod.config.tier.{tier}.{config.experiment}.simconfig.{simid}"
    sim_cfg = get_simconfig(config, tier, simid=simid)
    mac_subs = {}

    # configure generator
    if "generator" in sim_cfg:
        if not isinstance(sim_cfg.generator, str):
            msg = (
                "the field must be a string",
                f"{block}.generator",
            )
            raise SimflowConfigError(*msg)

        if sim_cfg.generator.startswith("~defines:"):
            key = sim_cfg.generator.removeprefix("~defines:")
            try:
                generator_lines = config.metadata.simprod.config.tier[tier][
                    config.experiment
                ].generators[key]
            except KeyError as e:
                msg = f"key {e} not found!"
                raise SimflowConfigError(msg, block) from e

        elif sim_cfg.generator.startswith("~vertices:"):
            vtx_file = patterns.vtx_filename_for_stp(config, simid, jobid="{JOBID}")
            generator_lines = [
                "/RMG/Generator/Confine FromFile",
                f"/RMG/Generator/Confinement/FromFile/FileName {nersc.dvs_ro(config, vtx_file)}",
            ]
            # in this case, vertex confinement is not required
            mac_subs["CONFINEMENT"] = None
        else:
            msg = (
                "the field must be prefixed with ~vertices: or ~defines:",
                f"{block}.generator",
            )
            raise SimflowConfigError(*msg)

        mac_subs["GENERATOR"] = "\n".join(generator_lines)

    # configure confinement
    if "confinement" in sim_cfg:
        confinement = None

        if isinstance(sim_cfg.confinement, str):
            if sim_cfg.confinement.startswith("~vertices:"):
                if sim_cfg.get("generator", "").startswith("~vertices:"):
                    msg = (
                        "no vertices in confinement field allowed if vertices are already specified as the generator",
                        f"{block}.confinement",
                    )
                    raise SimflowConfigError(*msg)

                vtx_file = patterns.vtx_filename_for_stp(config, simid, jobid="{JOBID}")
                confinement = [
                    "/RMG/Generator/Confine FromFile",
                    f"/RMG/Generator/Confinement/FromFile/FileName {nersc.dvs_ro(config, vtx_file)}",
                ]

            elif sim_cfg.confinement.startswith("~defines:"):
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
        try:
            text = lds.subst_vars(f.read().strip(), mac_subs, ignore_missing=False)
        except KeyError as e:
            msg = f"no rules found to substitute variable {e} in the macro template"
            raise SimflowConfigError(msg) from e

    # now write the macro to disk
    ofile = patterns.input_simjob_filename(config, tier=tier, simid=simid)
    ofile.parent.mkdir(parents=True, exist_ok=True)
    with ofile.open("w") as f:
        f.write(text)

    return text, ofile
