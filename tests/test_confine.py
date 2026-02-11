from __future__ import annotations

import pygeoml200
import pytest

from legendsimflow import commands, confine


@pytest.fixture(scope="session")
def test_generate_gdml(config):
    geom_config = config.metadata.simprod.config.geom["l200p02-geom-config"]

    return pygeoml200.core.construct(
        use_detailed_fiber_model=False, config=geom_config, public_geometry=True
    )


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
