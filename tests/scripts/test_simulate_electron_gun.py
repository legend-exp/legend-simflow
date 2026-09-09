from __future__ import annotations

import sys
from pathlib import Path

import lh5
import pytest
import yaml

from legendsimflow.metadata import ELECTRON_GUN_ENERGIES_IN_KEV
from legendsimflow.scripts import simulate_electron_gun

testprod = Path(__file__).parent.parent / "dummyprod"


@pytest.mark.needs_remage
def test_simulate_electron_gun(tmp_path, monkeypatch, l1000_config_factory):
    """The electron gun: geometry, macro, one remage run per energy, vertices plot."""
    config_path = l1000_config_factory(tmp_path)
    raw = yaml.safe_load(config_path.read_text())
    # the geometry generator and the macro template are read from the config
    # tree, which `$_` would otherwise resolve next to the temporary config
    raw["paths"]["config"] = str(testprod / "inputs/simprod/config")
    # run remage with few primaries
    raw["benchmark"] = {"enabled": True, "n_primaries": {"stp": 50}}
    config_path.write_text(yaml.safe_dump(raw))

    out_dir = tmp_path / "stp"
    geom_file = out_dir / "l1000dsg01-electron-gun-geom.gdml"
    macro_file = out_dir / "l1000dsg01-electron-gun-tier_stp.mac"
    stp_files = [
        out_dir / f"l1000dsg01-electron-gun-{e}keV-tier_stp.lh5"
        for e in ELECTRON_GUN_ENERGIES_IN_KEV[:2]
    ]
    plot_file = tmp_path / "plots" / "l1000dsg01-electron-gun-vertices.pdf"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "simulate-electron-gun",
            "--geom-config",
            str(testprod / "inputs/simprod/config/geom/l1000dsg01-geom-config.yaml"),
            "--geom-file",
            str(geom_file),
            "--macro-file",
            str(macro_file),
            "--stp-files",
            *[str(f) for f in stp_files],
            "--plot-file",
            str(plot_file),
            "--simflow-config",
            str(config_path),
        ],
    )
    simulate_electron_gun.main()

    assert geom_file.is_file()
    assert "/gps/particle e-" in macro_file.read_text()
    assert plot_file.is_file()
    for f in stp_files:
        assert f.is_file()
        assert "vtx" in lh5.ls(f)
        assert len(lh5.read("vtx", f)) == 50
        # electrons deposit energy in germanium only
        assert len(lh5.ls(f, "stp/")) > 0
