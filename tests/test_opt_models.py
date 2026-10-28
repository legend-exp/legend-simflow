from __future__ import annotations

import warnings

import numpy as np
import pyg4ometry
import pygeomtools
import pytest
from pyg4ometry.geant4 import LogicalVolume, PhysicalVolume, Registry, solid

from legendsimflow import opt_models


def _fiber_geometry(tmp_path, solids):
    reg = Registry()
    world = LogicalVolume(
        solid.Box("world", 10, 10, 10, reg, lunit="m"), "G4_Galactic", "world", reg
    )
    reg.setWorld(world)
    for i, make_solid in enumerate(solids):
        name = f"fiber_core_{i}"
        lv = LogicalVolume(make_solid(name, reg), "G4_POLYSTYRENE", name, reg)
        pv = PhysicalVolume([0, 0, 0], [i, 0, 0, "m"], lv, name, world, reg)
        pv.set_pygeom_active_detector(
            pygeomtools.RemageDetectorInfo(
                "scintillator", 1, allow_uid_reuse=True, ntuple_name="fiber"
            )
        )
    gdml_file = tmp_path / "geom.gdml"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pygeomtools.write_pygeom(reg, gdml_file)
    return pyg4ometry.gdml.Reader(gdml_file).getRegistry()


def _tubs(length):
    return lambda name, reg: solid.Tubs(name, 100, 101, length, 0, 1, reg, "mm")


def _box(length):
    return lambda name, reg: solid.Box(name, 1, 1, length, reg, "mm")


def test_uniform():
    x = np.zeros(3)
    assert opt_models.uniform(x, x, x, prob=0.02).tolist() == [0.02] * 3


def test_l1000_fibers(tmp_path):
    geometry = _fiber_geometry(tmp_path, [_tubs(1000), _box(1000)])
    # in the middle and at the top end of 1 m long fibers
    z_top = 349.092
    z = np.array([z_top - 500, z_top])
    pars = {"amplitudes": [0.2, 0.3, 0.5], "absorption_lengths_in_m": [0.05, 0.2, 4]}

    def end_prob(dist_in_mm):
        return sum(
            a * np.exp(-dist_in_mm / (1000 * length))
            for a, length in zip(*pars.values(), strict=True)
        )

    assert end_prob(0) == pytest.approx(1, rel=1e-3)
    prob = opt_models.l1000_fibers(z, z, z, **pars, geometry=geometry)
    assert prob == pytest.approx([2 * end_prob(500), end_prob(0) + end_prob(1000)])


@pytest.mark.parametrize(
    ("solids", "match"),
    [
        ([_tubs(1000), _box(900)], "single length"),
        (
            [_tubs(1000), lambda name, reg: solid.Orb(name, 1, reg, "mm")],
            "not a straight",
        ),
        ([], "single length"),
    ],
)
def test_l1000_fiber_length_errors(tmp_path, solids, match):
    with pytest.raises(ValueError, match=match):
        opt_models._l1000_fiber_length(_fiber_geometry(tmp_path, solids))
