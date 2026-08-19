from __future__ import annotations

import pytest

from legendsimflow import gamma_cascades


@pytest.fixture
def cascade_repo(tmp_path):
    """A minimal MAURINA-like repository with two isotopes."""
    for z, a, name in ((32, 76, "ge76"), (64, 155, "gd155")):
        d = tmp_path / str(z) / str(a)
        d.mkdir(parents=True)
        (d / f"{name}_ncapture_filelist.txt").write_text(
            f"cas{a}{name[:2]}_00100_a.txt 0 100000\n"
        )
    return tmp_path


@pytest.fixture
def cascade_config(fresh_config, cascade_repo):
    fresh_config.paths["maurina_gamma_cascades"] = str(cascade_repo)
    return fresh_config


def test_maurina_macro_commands(cascade_config, cascade_repo):
    lines = gamma_cascades.maurina_macro_commands(cascade_config)

    assert lines[0] == "/RMG/Processes/UseGrabmayrsGammaCascades true"
    assert lines[-1] == (
        "/RMG/GrabmayrGammaCascades/SetGammaCascadeRandomStartLocation 1"
    )

    # Z and A are taken from the two directory levels, sorted for determinism
    assert lines[1:-1] == [
        "/RMG/GrabmayrGammaCascades/SetGammaCascadeFilelist 32 76 "
        f"{cascade_repo / '32' / '76' / 'ge76_ncapture_filelist.txt'}",
        "/RMG/GrabmayrGammaCascades/SetGammaCascadeFilelist 64 155 "
        f"{cascade_repo / '64' / '155' / 'gd155_ncapture_filelist.txt'}",
    ]

    lines = gamma_cascades.maurina_macro_commands(cascade_config, subset={(32, 76)})

    assert lines[0] == "/RMG/Processes/UseGrabmayrsGammaCascades true"
    assert lines[-1] == (
        "/RMG/GrabmayrGammaCascades/SetGammaCascadeRandomStartLocation 1"
    )

    assert lines[1:-1] == [
        "/RMG/GrabmayrGammaCascades/SetGammaCascadeFilelist 32 76 "
        f"{cascade_repo / '32' / '76' / 'ge76_ncapture_filelist.txt'}",
    ]


def test_macro_block_empty_when_unconfigured(cascade_config):
    assert gamma_cascades.maurina_macro_block(cascade_config, {}) == ""
    assert (
        gamma_cascades.maurina_macro_block(
            cascade_config, {"maurina_gamma_cascades": False}
        )
        == ""
    )


def test_macro_block_registers_cascades(cascade_config):
    block = gamma_cascades.maurina_macro_block(
        cascade_config, {"maurina_gamma_cascades": True}
    )

    assert block.split("\n")[0] == "/RMG/Processes/UseGrabmayrsGammaCascades true"
