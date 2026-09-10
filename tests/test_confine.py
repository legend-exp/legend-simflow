from __future__ import annotations

import pyg4ometry
import pytest

from legendsimflow import commands, confine


def test_get_lar_minishroud_confine_commands(test_generate_gdml):
    lines = confine.get_lar_minishroud_confine_commands(test_generate_gdml, inside=True)

    assert len(lines) > 0
    assert isinstance(lines, list)

    assert "/RMG/Generator/Confinement/Geometrical/AddSolid Cylinder" in lines

    lines_outside = confine.get_lar_minishroud_confine_commands(
        test_generate_gdml, inside=False
    )

    assert len(lines_outside) > 0
    assert isinstance(lines_outside, list)

    assert (
        "/RMG/Generator/Confinement/Geometrical/AddExcludedSolid Cylinder"
        in lines_outside
    )

    with pytest.raises(ValueError):
        confine.get_lar_minishroud_confine_commands(
            test_generate_gdml, pattern="non_existent_pattern*"
        )
    with pytest.raises(ValueError):
        # exist pattern but not a nms
        confine.get_lar_minishroud_confine_commands(test_generate_gdml, pattern="V**")

    # test with eval

    lines_eval = commands.get_confinement_from_function(
        "legendsimflow.confine.get_lar_minishroud_confine_commands(<...>,inside=True)",
        test_generate_gdml,
    )
    assert lines_eval == lines

    lines_outside_eval = commands.get_confinement_from_function(
        "legendsimflow.confine.get_lar_minishroud_confine_commands(<...>,inside=False)",
        test_generate_gdml,
    )
    assert lines_outside_eval == lines_outside

    # test with string
    lines_eval_string = commands.get_confinement_from_function(
        "legendsimflow.confine.get_lar_minishroud_confine_commands(<...>,lar_name= 'liquid_argon',inside=True)",
        test_generate_gdml,
    )
    assert lines_eval_string == lines

    # without any kwarg
    lines_eval_args = commands.get_confinement_from_function(
        "legendsimflow.confine.get_lar_minishroud_confine_commands(<...>,'minishroud_tube*',True,lar_name= 'liquid_argon')",
        test_generate_gdml,
    )
    assert lines_eval_args == lines

    # outer cylinder: valid outside confinement
    lines_outer = confine.get_lar_minishroud_confine_commands(
        test_generate_gdml,
        inside=False,
        outer_radius_in_mm=1000,
        outer_height_in_mm=2000,
    )
    assert "/RMG/Generator/Confinement/Geometrical/AddSolid Cylinder" in lines_outer
    assert (
        "/RMG/Generator/Confinement/Geometrical/Cylinder/OuterRadius 1000 mm"
        in lines_outer
    )
    assert (
        "/RMG/Generator/Confinement/Geometrical/Cylinder/Height 2000 mm" in lines_outer
    )

    # outer cylinder: raises when inside=True
    with pytest.raises(ValueError):
        confine.get_lar_minishroud_confine_commands(
            test_generate_gdml,
            inside=True,
            outer_radius_in_mm=1000,
            outer_height_in_mm=2000,
        )

    # outer cylinder: raises when only one parameter is given
    with pytest.raises(ValueError):
        confine.get_lar_minishroud_confine_commands(
            test_generate_gdml, inside=False, outer_radius_in_mm=1000
        )
    with pytest.raises(ValueError):
        confine.get_lar_minishroud_confine_commands(
            test_generate_gdml, inside=False, outer_height_in_mm=2000
        )


def test_get_hpge_bulk_confine_commands(test_gdml_file):
    reg = pyg4ometry.gdml.Reader(str(test_gdml_file)).getRegistry()
    lines = confine.get_hpge_bulk_confine_commands(reg)

    assert lines[0] == "/RMG/Generator/Confine Volume"
    volumes = [
        line.removeprefix("/RMG/Generator/Confinement/Physical/AddVolume ")
        for line in lines[1:]
    ]
    assert len(volumes) > 0
    # one HPGe physical volume per line, sorted, no regex
    assert volumes == sorted(volumes)
    assert all(v[0] in "VPBC" and "*" not in v for v in volumes)

    lines_eval = commands.get_confinement_from_function(
        "legendsimflow.confine.get_hpge_bulk_confine_commands(<...>)", reg
    )
    assert lines_eval == lines
