from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

from . import SimflowConfig


def dvs_ro(
    config: SimflowConfig, path: str | Path | Sequence[str | Path]
) -> Path | list[Path]:
    """Turn ``/global/...`` file paths to ``/dvs_ro/...`` on NERSC.

    Note
    ----
    `config` must contain a ``nersc`` key mapped to a dictionary containing a
    ``dvs_ro: True`` key.
    """
    if isinstance(path, str):
        path = Path(path)

    if not config.nersc.dvs_ro:
        return path

    if isinstance(path, Path):
        return Path(re.sub("/global", "/dvs_ro", path.as_posix()))

    return [dvs_ro(config, p) for p in path]
