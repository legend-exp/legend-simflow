from __future__ import annotations

import shlex
from pathlib import Path

import pytest

from legendsimflow import SimflowConfigError, commands, patterns
from legendsimflow import metadata as mutils


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


def test_resolve_custom_g4ndl_spec():
    # plain string -> source only, no extra args
    spec = mutils.resolve_custom_g4ndl_spec("JEFF-3.3")
    assert spec["name"] == "JEFF-3.3"
    assert spec["source"] == "JEFF-3.3"
    assert spec["args"] == []

    # mapping with defaults only -> same name
    spec = mutils.resolve_custom_g4ndl_spec({"source": "JEFF-3.3"})
    assert spec["name"] == "JEFF-3.3"
    assert spec["args"] == []

    # options -> auto-disambiguated name + translated CLI args
    spec = mutils.resolve_custom_g4ndl_spec(
        {"source": "JEFF-3.3", "scale": 2.0, "no_substitution": True}
    )
    assert spec["name"] == "JEFF-3.3-scale2.0-nosubst"
    assert spec["args"] == ["--scale", "2.0", "--no-substitution"]

    # allow_incomplete is a readable token
    spec = mutils.resolve_custom_g4ndl_spec(
        {"source": "JEFF-3.3", "allow_incomplete": True}
    )
    assert spec["name"] == "JEFF-3.3-allowinc"

    # explicit name wins over auto-derivation
    spec = mutils.resolve_custom_g4ndl_spec(
        {"source": "JEFF-3.3", "scale": 2.0, "name": "mylib"}
    )
    assert spec["name"] == "mylib"

    # a custom substitution file cannot be auto-named -> requires explicit name
    with pytest.raises(SimflowConfigError):
        mutils.resolve_custom_g4ndl_spec(
            {"source": "JEFF-3.3", "substitution": "my_xs.dat"}
        )
    spec = mutils.resolve_custom_g4ndl_spec(
        {"source": "JEFF-3.3", "substitution": "my_xs.dat", "name": "jeff-ntof"}
    )
    assert spec["name"] == "jeff-ntof"
    assert spec["args"] == ["--substitution", "my_xs.dat"]

    # location-only options do not influence the name (they don't change the library)
    for neutral in ({"cache_dir": "/tmp/a"}, {"base_library": "/data/G4NDL4.7.1"}):
        assert (
            mutils.resolve_custom_g4ndl_spec({"source": "JEFF-3.3", **neutral})["name"]
            == "JEFF-3.3"
        )

    # errors: unknown option, missing source, path/URL source without a name
    with pytest.raises(SimflowConfigError):
        mutils.resolve_custom_g4ndl_spec({"source": "JEFF-3.3", "bogus": 1})
    with pytest.raises(SimflowConfigError):
        mutils.resolve_custom_g4ndl_spec({"scale": 2.0})
    with pytest.raises(SimflowConfigError):
        mutils.resolve_custom_g4ndl_spec({"source": "/data/G4NDL4.7"})
    assert (
        mutils.resolve_custom_g4ndl_spec(
            {"source": "/data/G4NDL4.7", "name": "g4ndl47"}
        )["name"]
        == "g4ndl47"
    )


def test_g4ndl_env_prefix():
    # bare binary: plain host variables only
    env = commands._g4ndl_env_prefix(["remage"], "/abs/lib")
    assert env[0] == "env"
    assert "G4NEUTRONHPDATA=/abs/lib" in env
    assert "G4PARTICLEHPDATA=/abs/lib" in env
    assert not any(e.startswith(("APPTAINERENV_", "SINGULARITYENV_")) for e in env)

    # apptainer: plain vars kept, plus APPTAINERENV_ forwarding (survives --cleanenv)
    env = commands._g4ndl_env_prefix(
        ["apptainer", "run", "--cleanenv", "remage.sif"], "/abs/lib"
    )
    assert "G4NEUTRONHPDATA=/abs/lib" in env
    assert "APPTAINERENV_G4NEUTRONHPDATA=/abs/lib" in env
    assert "APPTAINERENV_G4PARTICLEHPDATA=/abs/lib" in env

    # singularity: SINGULARITYENV_ forwarding
    env = commands._g4ndl_env_prefix(["singularity", "run", "remage.sif"], "/abs/lib")
    assert "SINGULARITYENV_G4NEUTRONHPDATA=/abs/lib" in env


def test_remage_cli_custom_g4ndl(fresh_config):
    config = fresh_config

    # without a custom_G4NDL field, remage is invoked directly (no env prefix)
    cmd = commands.remage_run(config, "birds_nest_K40", tier="stp")
    assert not cmd.startswith("env ")

    # opt in and check the env prefix points at the generated library
    config.metadata.simprod.config.tier.stp.legend.simconfig["birds_nest_K40"][
        "custom_G4NDL"
    ] = "JEFF-3.3"

    cmd = commands.remage_run(config, "birds_nest_K40", tier="stp")
    libdir = patterns.custom_g4ndl_dirname(config, "JEFF-3.3").resolve()
    tokens = shlex.split(cmd)
    assert tokens[0] == "env"
    assert f"G4NEUTRONHPDATA={libdir}" in tokens
    assert f"G4PARTICLEHPDATA={libdir}" in tokens
    assert "remage" in tokens
