from __future__ import annotations

import awkward as ak
import lh5
import numpy as np
import pytest

from legendsimflow import patterns
from legendsimflow.metadata import ELECTRON_GUN_ENERGIES_IN_KEV


@pytest.mark.needs_remage
def test_simulate_electron_gun(legend_electron_gun_paths):
    geom_file = legend_electron_gun_paths["geom_file"]
    stp_files = legend_electron_gun_paths["stp_files"]

    # the geometry is written out for the extraction step
    assert geom_file.exists()

    energies = [patterns.electron_gun_energy_from_path(f) for f in stp_files]
    assert energies == list(ELECTRON_GUN_ENERGIES_IN_KEV)

    for energy, stp_file in zip(energies, stp_files, strict=True):
        assert stp_file.exists()

        dets = {name.removeprefix("stp/") for name in lh5.ls(stp_file, "stp/")}
        # the electrons are confined in the bulk of every germanium volume
        assert dets >= {"V00001A", "V00001B", "V00002A", "V00002B"}

        steps = lh5.read("stp/V00001A", stp_file).view_as("ak")
        # electrons of a few MeV stop within a couple of mm: all the energy is
        # deposited in the detector they started in, but for the few events
        # losing part of it to escaping bremsstrahlung photons
        edep = ak.to_numpy(ak.sum(steps.edep, axis=-1))
        assert np.median(edep) == pytest.approx(energy, rel=1e-3)
        assert np.mean(edep > 0.99 * energy) > 0.8
