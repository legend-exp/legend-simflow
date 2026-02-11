from __future__ import annotations

import pyg4ometry
import pytest

from legendsimflow import commands, confine


def test_get_lar_minishroud_confine_commands(config):
    reg = pyg4ometry.gdml.Reader(config.paths.geom / "l200p02-geom.gdml").getRegistry()

    lines = confine.get_lar_minishroud_confine_commands(reg, inside=True)

    assert len(lines) > 0
    assert isinstance(lines, list)

    assert "/RMG/Generator/Confinement/Geometrical/AddSolid Cylinder" in lines

    lines_outside = confine.get_lar_minishroud_confine_commands(reg, inside=False)

    assert len(lines_outside) > 0
    assert isinstance(lines_outside, list)

    assert (
        "/RMG/Generator/Confinement/Geometrical/AddExcludedSolid Cylinder"
        in lines_outside
    )

    with pytest.raises(ValueError):
        confine.get_lar_minishroud_confine_commands(
            reg, pattern="non_existent_pattern*"
        )
    with pytest.raises(ValueError):
        # exist pattern but not a nms
        confine.get_lar_minishroud_confine_commands(reg, pattern="V**")

    # test with eval

    lines_eval = commands.get_confinement_from_function(
        "legendsimflow.confine.get_lar_minishroud_confine_commands(<...>,inside=True)",
        reg,
    )
    assert lines_eval == lines

    lines_outside_eval = commands.get_confinement_from_function(
        "legendsimflow.confine.get_lar_minishroud_confine_commands(<...>,inside=False)",
        reg,
    )
    assert lines_outside_eval == lines_outside

    lines_outside_eval = commands.get_confinement_from_function(
        "legendsimflow.confine.get_lar_minishroud_confine_commands(<...>,inside=False)",
        reg,
    )

    # test with string
    lines_eval_string = commands.get_confinement_from_function(
        "legendsimflow.confine.get_lar_minishroud_confine_commands(<...>,lar_name= 'liquid_argon',inside=True)",
        reg,
    )
    assert lines_eval_string == lines

    # without any kwarg
    lines_eval_args = commands.get_confinement_from_function(
        "legendsimflow.confine.get_lar_minishroud_confine_commands(<...>,'minishroud_tube*',True,lar_name= 'liquid_argon')",
        reg,
    )
    assert lines_eval_args == lines
