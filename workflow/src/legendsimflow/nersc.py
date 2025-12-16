from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path

from snakemake.io import InputFiles
from snakemake.script import Snakemake

from . import SimflowConfig


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
    # use the read-only path in all input items
    new_input_dict = {
        k: dvs_ro(snakemake.config, v) for k, v in snakemake.input.items()
    }

    # hack the inputs of the original snakemake object
    snakemake.input = InputFiles(fromdict=new_input_dict, plainstr=True)

    return snakemake
