from __future__ import annotations

import argparse
import re
import subprocess
from collections.abc import Iterable
from pathlib import Path

import yaml
from dbetto import AttrsDict
from legenddataflowscripts.workflow.utils import subst_vars
from legendmeta import LegendMetadata
from snakemake.io import InputFiles
from snakemake.script import Snakemake

from . import SimflowConfig, aggregate
from .exceptions import SimflowConfigError


def _partition(xs, n):
    k, r = divmod(len(xs), n)
    out, i = [], 0
    for j in range(n):
        s = k + (j < r)
        out.append(xs[i : i + s])
        i += s
    return out


def snakemake_nersc_cli():
    parser = argparse.ArgumentParser(
        description="Execute the Simflow on multiple nodes in parallel."
    )
    parser.add_argument(
        "-N", "--nodes", type=int, required=True, help="number of nodes"
    )
    args, extra = parser.parse_known_args()

    if args.nodes < 2:
        msg = "must parallelize over at least 2 nodes"
        raise ValueError(msg)

    cfg_path = Path("./simflow-config.yaml")
    if not cfg_path.is_file():
        msg = "this program must be executed in the directory where simflow-config.yaml resides"
        raise RuntimeError(msg)

    with cfg_path.open("r") as f:
        config = yaml.safe_load(f)

    subst_vars(
        config,
        var_values={"_": Path().resolve()},
        use_env=True,
        ignore_missing=False,
    )
    config = AttrsDict(config)

    # NOTE: this will attempt a clone of legend-metadata, if the directory does not exist
    metadata = LegendMetadata(config.paths.metadata, lazy=True)

    if "legend_metadata_version" in config:
        metadata.checkout(config.legend_metadata_version)

    config["metadata"] = metadata

    simlist = config.get("simlist", None)
    make_tiers = config.make_tiers
    if simlist is None:
        # auto determine tier from config
        tiers = ("pdf", "cvt", "evt", "hit", "opt", "stp")
        tier = next(t for t in tiers if t in make_tiers)

        simlist = [
            f"{tier}.{simid}" for simid in aggregate.gen_list_of_all_simids(config)
        ]

    procs = []
    for simlist_chunk in _partition(simlist, args.nodes):
        smk_cmd = [
            "srun",
            "--nodes",
            "1",
            "--ntasks",
            "1",
            "--cpus-per-task",
            "256",
            "snakemake",
            "--worflow-profile",
            "workflow/profiles/nersc",
            "--config",
            "simlist=" + ",".join(simlist_chunk),
            *extra,
        ]

        print("INFO: spawning process:", " ".join(smk_cmd))  # noqa: T201
        procs.append(subprocess.Popen(smk_cmd))

    for p in procs:
        rc = p.wait()
        if rc != 0:
            msg = f"process failed: {p.args}"
            raise RuntimeError(msg)

    print("INFO: all snakemake processes successfully returned")  # noqa: T201


def dvs_ro(
    config: SimflowConfig, path: str | Path | Iterable[str | Path]
) -> str | Path | list[str | Path]:
    """Turn ``/global/...`` file paths to ``/dvs_ro/...`` on NERSC.

    The input type is preserved.

    Note
    ----
    `config` must contain a ``nersc`` key mapped to a dictionary containing a
    ``dvs_ro: True`` key.
    """
    if not config.nersc.dvs_ro:
        return path

    if isinstance(path, str):
        return re.sub("/global", "/dvs_ro", path)

    if isinstance(path, Path):
        return Path(re.sub("/global", "/dvs_ro", path.as_posix()))

    return [dvs_ro(config, p) for p in path]


def dvs_ro_snakemake(snakemake: Snakemake) -> Snakemake:
    """Swap the read-only filesystem path in all Snakemake input files.

    This function is meant to be used in Snakemake scripts, where the Snakemake
    rule attributes (input, output, ...) are accessible from the special object
    ``snakemake``.

    Warning
    -------
    This function mutates the input `snakemake` object in place.

    See also
    --------
    dvs_ro
    """
    # is this a named list?
    if len(snakemake.input.keys()) != 0:
        # use the read-only path in all input items
        new_input_dict = {
            k: dvs_ro(snakemake.config, v) for k, v in snakemake.input.items()
        }

        # hack the inputs of the original snakemake object
        snakemake.input = InputFiles(fromdict=new_input_dict, plainstr=True)
    # or an un-named list?
    else:
        new_input_list = [dvs_ro(snakemake.config, v) for v in snakemake.input]
        snakemake.input = InputFiles(toclone=new_input_list, plainstr=True)

    return snakemake


def is_scratch_enabled(config: SimflowConfig) -> bool:
    """Is the scratch folder enabled in this workflow?"""
    field = config.nersc.scratch

    if isinstance(field, bool) and not field:
        return False

    if not isinstance(field, str):
        msg = (
            "this field can be either false or path to a scratch folder",
            "config.nersc.scratch",
        )
        raise SimflowConfigError(*msg)

    if Path(field).exists() and not Path(field).is_dir():
        msg = (f"{field}: not a directory", "config.nersc.scratch")
        raise SimflowConfigError(*msg)

    return True


def scratch_dir(config: SimflowConfig) -> Path:
    """The scratch folder path configured in this workflow."""
    if not is_scratch_enabled(config):
        msg = "scratch folder not set"
        raise RuntimeError(msg)

    return Path(config.nersc.scratch)


def on_scratch(config: SimflowConfig, path: str | Path) -> Path:
    """Return the path of the file in the scratch folder."""
    path = Path(path)

    # rebuild path from parts to normalize (removes //, ./, etc.)
    parts = path.parts

    # drop anchor if absolute
    if path.is_absolute():
        parts = parts[1:]

    path = Path(*parts)

    new_path = scratch_dir(config) / path
    new_path.parent.mkdir(parents=True, exist_ok=True)
    return new_path
