from __future__ import annotations

import re
import shutil
from collections.abc import Callable, Iterable
from pathlib import Path

from snakemake.io import InputFiles
from snakemake.script import Snakemake

from . import SimflowConfig
from .exceptions import SimflowConfigError


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

    return Path(config.nersc.scratch) / config._proctime


def on_scratch(config: SimflowConfig, path: str | Path) -> Path:
    """Return the path of the file in the scratch folder.

    Also makes sure the parent folder exists.
    """
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


def make_on_scratch(config: SimflowConfig, path: str | Path) -> tuple[Path, Callable]:
    """Return tools to produce a file in the scratch folder and move it to the final destination.

    Returns the temporary path in the scratch dir and a function that will move
    it back to the original destination.
    """
    if not is_scratch_enabled(config):
        return path, lambda: None

    scratch_path = on_scratch(config, path)

    def move(scratch_path=scratch_path, path=path):
        shutil.move(scratch_path, path)

    return scratch_path, move
