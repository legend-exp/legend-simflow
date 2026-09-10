from __future__ import annotations

import shlex
from pathlib import Path

import pytest

from legendsimflow import SimflowConfigError, commands, patterns


def test_confine_by_volume_commands():
    # Test with surface sampling enabled
    lines = commands._confine_by_volume(
        is_surface=True, volume="test_volume", surface_max_intersections=50
    )
    assert "/RMG/Generator/Confinement/Physical/AddVolume test_volume" in lines
    assert "/RMG/Generator/Confinement/SampleOnSurface true" in lines
    assert "/RMG/Generator/Confinement/SurfaceSampleMaxIntersections 50" in lines

    # Test with surface sampling disabled
    lines = commands._confine_by_volume(is_surface=False, volume="test_volume")
    assert "/RMG/Generator/Confinement/Physical/AddVolume test_volume" in lines
    assert "/RMG/Generator/Confinement/SampleOnSurface true" not in lines
    assert "/RMG/Generator/Confinement/SurfaceSampleMaxIntersections 100" not in lines


def test_make_macro(config):
    text, fmac = commands.make_remage_macro(config, "birds_nest_K40", "stp")

    assert (
        fmac
        == Path(config.paths.macros)
        / f"{config.experiment}-birds_nest_K40-tier_stp.mac"
    )
    assert fmac.is_file()

    assert set(
        config.metadata.simprod.config.tier.stp.legend.confinement["birds_nest"]
    ).issubset(text.split("\n"))
    assert set(
        config.metadata.simprod.config.tier.stp.legend.generators["K40"]
    ).issubset(text.split("\n"))

    text, fmac = commands.make_remage_macro(config, "pen_plates_Ra224_to_Pb208", "stp")
    assert set(
        config.metadata.simprod.config.tier.stp.legend.generators["Ra224_to_Pb208"]
    ).issubset(text.split("\n"))

    confine = [
        "/RMG/Generator/Confine Volume",
        "/RMG/Generator/Confinement/Physical/AddVolume hpge_assembly_plate_pen.*",
    ]
    assert set(confine).issubset(text.split("\n"))
    assert "/RMG/Generator/Confinement/SampleOnSurface" not in text

    text, fmac = commands.make_remage_macro(
        config, "phbr_surface_Ra228_to_Ac228", "stp"
    )
    confine = [
        "/RMG/Generator/Confine Volume",
        "/RMG/Generator/Confinement/Physical/AddVolume hpge_assembly_phbr_spring.*",
        "/RMG/Generator/Confinement/Physical/AddVolume hpge_assembly_phbr_washer.*",
        "/RMG/Generator/Confinement/SampleOnSurface true",
        "/RMG/Generator/Confinement/SurfaceSampleMaxIntersections 100",
    ]
    assert set(confine).issubset(text.split("\n"))

    text, fmac = commands.make_remage_macro(
        config, "hpge_bulk_high_thr_Rn222_to_Po214", "stp"
    )

    assert text is not None

    text, fmac = commands.make_remage_macro(config, "lar_hpge_shell_K42", "stp")
    confine = [
        "/RMG/Generator/Confine FromFile",
        "/RMG/Generator/Confinement/FromFile/FileName "
        + str(
            patterns.vtx_filename_for_stp(config, "lar_hpge_shell_K42", jobid="{JOBID}")
        ),
    ]
    assert set(confine).issubset(text.split("\n"))

    text, fmac = commands.make_remage_macro(config, "exotic_physics_process", "stp")
    confine = [
        "/RMG/Generator/Confine FromFile",
        "/RMG/Generator/Confinement/FromFile/FileName "
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
        "/RMG/Generator/Confinement/FromFile/FileName "
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

    metadata.simprod.config.tier.stp.legend.simconfig["birds_nest_K40"]["generator"] = (
        "coddue"
    )
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")

    metadata.simprod.config.tier.stp.legend.simconfig["birds_nest_K40"]["generator"] = (
        "~coddue:boh"
    )
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")

    metadata.simprod.config.tier.stp.legend.simconfig["birds_nest_K40"]["generator"] = (
        "~defines:boh"
    )
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")


def test_make_macro_errors_2(fresh_config):
    config = fresh_config
    metadata = fresh_config.metadata

    metadata.simprod.config.tier.stp.legend.simconfig["birds_nest_K40"][
        "confinement"
    ] = "~baaaaaa:beh"

    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")

    metadata.simprod.config.tier.stp.legend.simconfig["birds_nest_K40"][
        "confinement"
    ] = "~defines:beh"
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")

    metadata.simprod.config.tier.stp.legend.simconfig["birds_nest_K40"][
        "confinement"
    ] = {}
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")


def test_make_macro_errors_confinement_list(fresh_config):
    """An invalid entry anywhere in a confinement list gives a config error."""
    config = fresh_config
    metadata = fresh_config.metadata
    simconfig = metadata.simprod.config.tier.stp.legend.simconfig["birds_nest_K40"]

    # invalid entry first: the valid one that follows must not resurrect the
    # accumulator and crash with a TypeError instead
    simconfig["confinement"] = [
        "~baaaaaa:beh",
        "~volumes.bulk:hpge_assembly_plate_pen.*",
    ]
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")

    # invalid entry last
    simconfig["confinement"] = [
        "~volumes.bulk:hpge_assembly_plate_pen.*",
        "~baaaaaa:beh",
    ]
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "birds_nest_K40", "stp")


def test_make_macro_errors_vertices(fresh_config):
    config = fresh_config
    metadata = fresh_config.metadata

    metadata.simprod.config.tier.stp.legend.simconfig.exotic_physics_process[
        "confinement"
    ] = "~vertices:blah"
    with pytest.raises(SimflowConfigError):
        commands.make_remage_macro(config, "exotic_physics_process", "stp")

    metadata.simprod.config.tier.stp.legend.simconfig.exotic_physics_process.pop(
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

    cmd = commands.remage_run(
        config, "birds_nest_K40", tier="stp", geom="/some/geom.gdml", macro_free=True
    )
    mac_cmds = shlex.split(cmd.partition(" -- ")[2])
    assert all(cmd[0] == "/" for cmd in mac_cmds)

    # macro_free renders the macro here and now, so geom must be a real path
    with pytest.raises(ValueError, match="Snakemake placeholder"):
        commands.remage_run(config, "birds_nest_K40", tier="stp", macro_free=True)

    # a path that merely contains braces is not a placeholder
    cmd = commands.remage_run(
        config,
        "birds_nest_K40",
        tier="stp",
        geom="/some/geom{v1}.gdml",
        macro_free=True,
    )
    assert "/some/geom{v1}.gdml" in shlex.split(cmd.partition(" -- ")[0])

    config.benchmark.enabled = True
    config.benchmark.n_primaries.stp = 999

    cmd = commands.remage_run(config, "birds_nest_K40", tier="stp")
    cmdline = shlex.split(cmd.partition(" -- ")[0])
    assert "N_EVENTS=999" in cmdline

    cmd = commands.remage_run(config, "lar_hpge_shell_K42", tier="stp", jobid="0001")
    cmdline = shlex.split(cmd.partition(" -- ")[0])
    assert "JOBID=0001" in cmdline


def test_remage_cli_scratch_mv_is_quoted(fresh_config, tmp_path):
    """The scratch-to-final move must survive paths with shell metacharacters."""
    config = fresh_config
    scratch = tmp_path / "scratch dir"
    config.nersc.scratch = str(scratch)

    cmd = commands.remage_run(config, "birds_nest_K40", tier="stp")

    remage_part, _, mv_part = cmd.partition(" && ")
    assert mv_part.startswith("mv -v ")

    # the whole command must still tokenize, with the two paths intact
    mv_tokens = shlex.split(mv_part)
    assert mv_tokens[:2] == ["mv", "-v"]
    assert len(mv_tokens) == 4

    src, dest = Path(mv_tokens[2]), Path(mv_tokens[3])
    assert src.is_relative_to(scratch)
    assert not dest.is_relative_to(scratch)
    # the remage invocation writes to the scratch copy
    assert str(src) in shlex.split(remage_part)
