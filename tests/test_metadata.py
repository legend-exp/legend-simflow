from __future__ import annotations

from dbetto import AttrsDict

from legendsimflow import metadata


def test_all(config):
    assert isinstance(
        metadata.get_vtx_simconfig(config, "lar_hpge_shell_K42"), AttrsDict
    )
