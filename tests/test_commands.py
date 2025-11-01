from __future__ import annotations

import shlex
from pathlib import Path

import pytest

from legendsimflow import SimflowConfigError, commands, patterns


def test_make_macro(config):
    text, fmac = commands.make_remage_macro(config, "birds-nest-K40", "stp")

    assert (
        fmac
        == Path(config.paths.macros)
        / f"stp/{config.experiment}-birds-nest-K40-tier_stp.mac"
    )
    assert fmac.is_file()

    assert set(
        config.metadata.simprod.config.tier.stp.l200p03.confinement["birds-nest"]
    ).issubset(text.split("\n"))
    assert set(
        config.metadata.simprod.config.tier.stp.l200p03.generators["K40"]
    ).issubset(text.split("\n"))

    text, fmac = commands.make_remage_macro(config, "pen-plates-Ra224-to-Pb208", "stp")
    assert set(
        config.metadata.simprod.config.tier.stp.l200p03.generators["Ra224-to-Pb208"]
    ).issubset(text.split("\n"))

    confine = [
        "/RMG/Generator/Confine Volume",
        "/RMG/Generator/Confinement/Physical/AddVolume pen.*",
    ]
    assert set(confine).issubset(text.split("\n"))
    assert "/RMG/Generator/Confinement/SampleOnSurface" not in text

    text, fmac = commands.make_remage_macro(
        config, "phbr-surface-Ra228-to-Ac228", "stp"
    )
    confine = [
        "/RMG/Generator/Confine Volume",
        "/RMG/Generator/Confinement/Physical/AddVolume phbr_spring.*",
        "/RMG/Generator/Confinement/Physical/AddVolume phbr_washer.*",
        "/RMG/Generator/Confinement/SampleOnSurface true",
    ]
    assert set(confine).issubset(text.split("\n"))

    text, fmac = commands.make_remage_macro(
        config, "hpge-bulk-high-thr-Rn222-to-Po214", "stp"
    )
    assert text is not None


def test_make_macro_errors_1(fresh_config):
    config = fresh_config
    metadata = fresh_config.metadata

    metadata.simprod.config.tier.stp.l200p03.simconfig["birds-nest-K40"][
        "generator"
    ] = "coddue"
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds-nest-K40", "stp")

    metadata.simprod.config.tier.stp.l200p03.simconfig["birds-nest-K40"][
        "generator"
    ] = "~coddue:boh"
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds-nest-K40", "stp")

    metadata.simprod.config.tier.stp.l200p03.simconfig["birds-nest-K40"][
        "generator"
    ] = "~defines:boh"
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds-nest-K40", "stp")


def test_make_macro_errors_2(fresh_config):
    config = fresh_config
    metadata = fresh_config.metadata

    metadata.simprod.config.tier.stp.l200p03.simconfig["birds-nest-K40"][
        "confinement"
    ] = "~baaaaaa:beh"

    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds-nest-K40", "stp")

    metadata.simprod.config.tier.stp.l200p03.simconfig["birds-nest-K40"][
        "confinement"
    ] = "~defines:beh"
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds-nest-K40", "stp")

    metadata.simprod.config.tier.stp.l200p03.simconfig["birds-nest-K40"][
        "confinement"
    ] = {}
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds-nest-K40", "stp")


def test_remage_cli(fresh_config):
    config = fresh_config

    cmd = commands.remage_run(config, "birds-nest-K40", "stp")
    assert isinstance(cmd, str)
    assert len(cmd) > 0
    assert shlex.split(cmd)[-1] == patterns.input_simjob_filename(
        config, tier="stp", simid="birds-nest-K40"
    )

    cmd = commands.remage_run(config, "birds-nest-K40", "stp", macro_free=True)
    mac_cmds = shlex.split(cmd.partition(" -- ")[2])
    assert all(cmd[0] == "/" for cmd in mac_cmds)

    config.benchmark.enabled = True
    config.benchmark.n_primaries.stp = 999

    cmd = commands.remage_run(config, "birds-nest-K40", "stp")
    cmdline = shlex.split(cmd.partition(" -- ")[0])
    assert "N_EVENTS=999" in cmdline
