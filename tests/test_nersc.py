from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from legendsimflow import SimflowConfigError, nersc


def test_dvs_ro(fresh_config):
    c = fresh_config
    c.nersc.dvs_ro = True

    assert nersc.dvs_ro(c, "/global/test/file") == "/dvs_ro/test/file"
    assert nersc.dvs_ro(c, Path("/global/test/file")) == Path("/dvs_ro/test/file")
    assert nersc.dvs_ro(c, "/another/test/file") == "/another/test/file"
    assert nersc.dvs_ro(c, ["/global/test/file", "/another/test/file"]) == [
        "/dvs_ro/test/file",
        "/another/test/file",
    ]

    c.nersc.dvs_ro = False
    assert nersc.dvs_ro(c, "/global/test/file") == "/global/test/file"


def test_scratch(fresh_config):
    c = fresh_config

    c.nersc.scratch = "/tmp"
    assert nersc.is_scratch_enabled(c) is True

    c.nersc.scratch = False
    assert nersc.is_scratch_enabled(c) is False

    c.nersc.scratch = True
    with pytest.raises(SimflowConfigError):
        nersc.is_scratch_enabled(c)

    c.nersc.scratch = __file__
    with pytest.raises(SimflowConfigError):
        nersc.is_scratch_enabled(c)

    c.nersc.scratch = False
    with pytest.raises(RuntimeError):
        nersc.scratch_dir(c)

    tmp = Path(tempfile.gettempdir())
    c.nersc.scratch = str(tmp)
    assert nersc.scratch_dir(c) == tmp

    c.nersc.scratch = str(tmp / "legend-simflow")
    assert (
        nersc.on_scratch(c, "test.ext").as_posix() == f"{tmp}/legend-simflow/test.ext"
    )
    assert Path(f"{tmp}/legend-simflow").is_dir()

    assert nersc.on_scratch(c, __file__).as_posix() == f"{tmp}/legend-simflow{__file__}"
