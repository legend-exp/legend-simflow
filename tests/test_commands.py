from __future__ import annotations

import shlex
from pathlib import Path

import pytest

from legendsimflow import SimflowConfigError, commands, patterns


def test_make_macro(config):
    text, fmac = commands.make_remage_macro(config, "birds_nest_K40", "stp")

    assert (
        fmac
        == Path(config.paths.macros)
        / f"stp/{config.experiment}-birds_nest_K40-tier_stp.mac"
    )
    assert fmac.is_file()

    assert set(
        config.metadata.simprod.config.tier.stp.l200p02.confinement["birds_nest"]
    ).issubset(text.split("\n"))
    assert set(
        config.metadata.simprod.config.tier.stp.l200p02.generators["K40"]
    ).issubset(text.split("\n"))

    text, fmac = commands.make_remage_macro(config, "pen_plates_Ra224_to_Pb208", "stp")
    assert set(
        config.metadata.simprod.config.tier.stp.l200p02.generators["Ra224_to_Pb208"]
    ).issubset(text.split("\n"))

    confine = [
        "/RMG/Generator/Confine Volume",
        "/RMG/Generator/Confinement/Physical/AddVolume pen.*",
    ]
    assert set(confine).issubset(text.split("\n"))
    assert "/RMG/Generator/Confinement/SampleOnSurface" not in text

    text, fmac = commands.make_remage_macro(
        config, "phbr_surface_Ra228_to_Ac228", "stp"
    )
    confine = [
        "/RMG/Generator/Confine Volume",
        "/RMG/Generator/Confinement/Physical/AddVolume phbr_spring.*",
        "/RMG/Generator/Confinement/Physical/AddVolume phbr_washer.*",
        "/RMG/Generator/Confinement/SampleOnSurface true",
    ]
    assert set(confine).issubset(text.split("\n"))

    text, fmac = commands.make_remage_macro(
        config, "hpge_bulk_high_thr_Rn222_to_Po214", "stp"
    )
    assert text is not None

    text, fmac = commands.make_remage_macro(config, "lar_hpge_shell_K42", "stp")
    confine = [
        "/RMG/Generator/Confine FromFile",
        "/RMG/Generator/FromFile/FileName "
        + str(
            patterns.vtx_filename_for_stp(config, "lar_hpge_shell_K42", jobid="{JOBID}")
        ),
    ]
    assert set(confine).issubset(text.split("\n"))

    text, fmac = commands.make_remage_macro(config, "exotic_physics_process", "stp")
    confine = [
        "/RMG/Generator/Confine FromFile",
        "/RMG/Generator/FromFile/FileName "
        + str(
            patterns.vtx_filename_for_stp(
                config, "exotic_physics_process", jobid="{JOBID}"
            )
        ),
    ]
    assert set(confine).issubset(text.split("\n"))
    assert "/RMG/Generator/Select" not in text

    text, fmac = commands.make_remage_macro(config, "exotic_physics_hpge", "stp")
    confine = [
        "/RMG/Generator/Confine FromFile",
        "/RMG/Generator/FromFile/FileName "
        + str(
            patterns.vtx_filename_for_stp(
                config, "exotic_physics_hpge", jobid="{JOBID}"
            )
        ),
    ]
    assert set(confine).issubset(text.split("\n"))
    assert "/RMG/Generator/Confine Volume" in text


def test_make_macro_errors_1(fresh_config):
    config = fresh_config
    metadata = fresh_config.metadata

    metadata.simprod.config.tier.stp.l200p02.simconfig["birds_nest_K40"][
        "generator"
    ] = "coddue"
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")

    metadata.simprod.config.tier.stp.l200p02.simconfig["birds_nest_K40"][
        "generator"
    ] = "~coddue:boh"
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")

    metadata.simprod.config.tier.stp.l200p02.simconfig["birds_nest_K40"][
        "generator"
    ] = "~defines:boh"
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")


def test_make_macro_errors_2(fresh_config):
    config = fresh_config
    metadata = fresh_config.metadata

    metadata.simprod.config.tier.stp.l200p02.simconfig["birds_nest_K40"][
        "confinement"
    ] = "~baaaaaa:beh"

    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")

    metadata.simprod.config.tier.stp.l200p02.simconfig["birds_nest_K40"][
        "confinement"
    ] = "~defines:beh"
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")

    metadata.simprod.config.tier.stp.l200p02.simconfig["birds_nest_K40"][
        "confinement"
    ] = {}
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")


def test_make_macro_errors_vertices(fresh_config):
    config = fresh_config
    metadata = fresh_config.metadata

    metadata.simprod.config.tier.stp.l200p02.simconfig.exotic_physics_process[
        "confinement"
    ] = "~vertices:blah"
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "exotic_physics_process", "stp")

    metadata.simprod.config.tier.stp.l200p02.simconfig.exotic_physics_process.pop(
        "generator"
    )
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "exotic_physics_process", "stp")


def test_remage_cli(fresh_config):
    config = fresh_config

    cmd = commands.remage_run(config, "birds_nest_K40", tier="stp")
    assert isinstance(cmd, str)
    assert len(cmd) > 0
    assert (
        shlex.split(cmd)[-1]
        == patterns.input_simjob_filename(
            config, tier="stp", simid="birds_nest_K40"
        ).as_posix()
    )

    cmd = commands.remage_run(config, "birds_nest_K40", tier="stp", macro_free=True)
    mac_cmds = shlex.split(cmd.partition(" -- ")[2])
    assert all(cmd[0] == "/" for cmd in mac_cmds)

    config.benchmark.enabled = True
    config.benchmark.n_primaries.stp = 999

    cmd = commands.remage_run(config, "birds_nest_K40", tier="stp")
    cmdline = shlex.split(cmd.partition(" -- ")[0])
    assert "N_EVENTS=999" in cmdline

    cmd = commands.remage_run(config, "lar_hpge_shell_K42", tier="stp", jobid="0001")
    cmdline = shlex.split(cmd.partition(" -- ")[0])
    assert "JOBID=0001" in cmdline
