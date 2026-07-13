from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml
from lgdo import Array, Table, lh5

from legendsimflow.scripts import get_aoe_energy_correction

dummyprod = Path(__file__).parent.parent / "dummyprod"


def _write_fake_mc_hit(path: Path, det: str = "V05261A", n_events: int = 8000) -> None:
    rng = np.random.default_rng(12345)
    energy = rng.uniform(900.0, 2400.0, n_events)

    def _make_aoe(sigma: float, frac_bkg: float = 0.15) -> np.ndarray:
        # realistic A/E: a signal Gaussian around 1.0 plus a low-side tail
        # (multi-site background), as expected by pygama's compton-band fitter
        signal = rng.normal(1.0, sigma, n_events)
        tail = 1.0 - rng.exponential(4.0 * sigma, n_events)
        is_bkg = rng.random(n_events) < frac_bkg
        return np.where(is_bkg, tail, signal)

    mc_hit = Table(
        col_dict={
            "psd": Table(
                col_dict={
                    "single_temp": Table(col_dict={"aoe_raw": Array(_make_aoe(0.01))}),
                    "pulse_lib": Table(col_dict={"aoe_raw": Array(_make_aoe(0.012))}),
                }
            ),
            "energy": Array(energy),
        }
    )
    lh5.write(mc_hit, f"hit/{det}", path, wo_mode="write_safe")


def _build_argv(tmp_path: Path) -> list[str]:
    config_path = tmp_path / "simflow-config.yaml"
    raw = yaml.safe_load((dummyprod / "simflow-config.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    raw["paths"]["l200data"] = str(tmp_path / "missing-l200data")
    raw["paths"]["tier_hit"] = str(tmp_path)
    config_path.write_text(yaml.safe_dump(raw))

    return [
        "get-aoe-energy-correction",
        "--hpge-detector",
        "V05261A",
        "--pars-file",
        str(tmp_path / "outputs" / "aoe-energy-correction.yaml"),
        "--plot-file",
        str(tmp_path / "plots" / "aoe-energy-correction.pdf"),
        "--simflow-config",
        str(config_path),
    ]


def test_get_aoe_energy_correction_mc_only(tmp_path, monkeypatch):
    simid_dir = tmp_path / "generated" / "tier" / "hit" / "test-Pb212"
    simid_dir.mkdir(parents=True, exist_ok=True)
    _write_fake_mc_hit(simid_dir / "test-hit.lh5")

    monkeypatch.setattr(sys, "argv", _build_argv(tmp_path))

    get_aoe_energy_correction.main()

    pars_file = tmp_path / "outputs" / "aoe-energy-correction.yaml"
    plot_file = tmp_path / "plots" / "aoe-energy-correction.pdf"
    assert pars_file.exists()
    assert plot_file.exists()

    out = yaml.safe_load(pars_file.read_text())
    assert set(out) == {"energy_corrections", "mc_resolution"}
    assert set(out["energy_corrections"]) == {"single_template", "psl"}
    assert "data_validation" not in out

    # the mean energy correction is the primary product: its fit parameters must
    # be real numbers (a NaN would mean the underlying fit silently failed)
    for key in ("single_template", "psl"):
        pars = out["energy_corrections"][key]["pars"]
        assert np.isfinite(pars["a"])
        assert np.isfinite(pars["b"])
