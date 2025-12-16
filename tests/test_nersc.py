from __future__ import annotations

from pathlib import Path

from legendsimflow import nersc


def test_nersc_utils(fresh_config):
    c = fresh_config
    c.nersc.dvs_ro = True

    assert nersc.dvs_ro(c, "/global/test/file") == Path("/dvs_ro/test/file")
    assert nersc.dvs_ro(c, Path("/global/test/file")) == Path("/dvs_ro/test/file")
    assert nersc.dvs_ro(c, "/another/test/file") == Path("/another/test/file")
    assert nersc.dvs_ro(c, ["/global/test/file", "/another/test/file"]) == [
        Path("/dvs_ro/test/file"),
        Path("/another/test/file"),
    ]

    c.nersc.dvs_ro = False
    assert nersc.dvs_ro(c, "/global/test/file") == Path("/global/test/file")
