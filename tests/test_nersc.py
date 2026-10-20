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

    # must not replace /global not at the start of the path
    assert nersc.dvs_ro(c, "/some/global/test/file") == "/some/global/test/file"
    assert nersc.dvs_ro(c, Path("/some/global/test/file")) == Path(
        "/some/global/test/file"
    )

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
    assert nersc.scratch_dir(c).parent == tmp

    c.nersc.scratch = str(tmp)
    scr_path = nersc.on_scratch(c, "test.ext")
    assert scr_path.name == "test.ext"
    assert scr_path.parent.parent == Path(c.nersc.scratch)

    scr_path = nersc.on_scratch(c, "/some/test.ext", run_folder=False)
    assert scr_path == tmp / "some/test.ext"


def test_make_on_scratch_copies(fresh_config, tmp_path):
    # a move within one filesystem is a rename and keeps the inode, and with it
    # the group the file had on scratch, defeating the setgid bit of the
    # destination folder
    c = fresh_config
    c.nersc.scratch = str(tmp_path / "scratch")

    final = tmp_path / "out/file.ext"
    final.parent.mkdir(parents=True)

    scratch_path, copy2final = nersc.make_on_scratch(c, final)
    assert scratch_path != final
    scratch_path.write_text("data")
    scratch_inode = scratch_path.stat().st_ino

    copy2final()
    assert final.read_text() == "data"
    assert final.stat().st_ino != scratch_inode
    assert not scratch_path.exists()


def test_make_on_scratch_without_scratch(fresh_config, tmp_path):
    c = fresh_config
    c.nersc.scratch = False

    final = tmp_path / "file.ext"
    scratch_path, copy2final = nersc.make_on_scratch(c, final)
    assert scratch_path == final

    final.write_text("data")
    copy2final()
    assert final.read_text() == "data"
