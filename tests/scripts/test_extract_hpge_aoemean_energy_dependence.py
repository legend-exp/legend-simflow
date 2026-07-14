from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path

import lh5
import numpy as np
import pytest
import yaml
from lgdo import Array, Table

from legendsimflow.scripts import extract_hpge_aoemean_energy_dependence

dummyprod = Path(__file__).parent.parent / "dummyprod"

# per-method A/E width used to synthesize the fake hit data
_SIM_TYPE_SIGMA = {"single_temp": 0.01, "pulse_lib": 0.012}


def _write_fake_mc_hit(
    path: Path,
    det: str = "V05261A",
    n_events: int = 40000,
    sim_types: Sequence[str] = ("single_temp", "pulse_lib"),
) -> None:
    rng = np.random.default_rng(12345)
    energy = rng.uniform(900.0, 2400.0, n_events)

    def _make_aoe(sigma: float, frac_bkg: float = 0.15) -> np.ndarray:
        # realistic A/E: a signal Gaussian around 1.0 plus a low-side tail
        # (multi-site background), as expected by pygama's compton-band fitter
        signal = rng.normal(1.0, sigma, n_events)
        tail = 1.0 - rng.exponential(4.0 * sigma, n_events)
        is_bkg = rng.random(n_events) < frac_bkg
        return np.where(is_bkg, tail, signal)

    psd = Table(
        col_dict={
            st: Table(col_dict={"aoe_raw": Array(_make_aoe(_SIM_TYPE_SIGMA[st]))})
            for st in sim_types
        }
    )
    mc_hit = Table(col_dict={"psd": psd, "energy": Array(energy)})
    lh5.write(mc_hit, f"hit/{det}", path, wo_mode="write_safe")


def _build_argv(tmp_path: Path, hit_files: list[Path]) -> list[str]:
    config_path = tmp_path / "simflow-config.yaml"
    raw = yaml.safe_load((dummyprod / "simflow-config.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    raw["paths"]["l200data"] = str(tmp_path / "missing-l200data")
    config_path.write_text(yaml.safe_dump(raw))

    return [
        "extract-hpge-aoemean-energy-dependence",
        "--hpge-detector",
        "V05261A",
        "--hit-files",
        *[str(f) for f in hit_files],
        "--pars-file",
        str(tmp_path / "outputs" / "V05261A-model.yaml"),
        "--plot-file",
        str(tmp_path / "plots" / "V05261A-fit-result.pdf"),
        "--simflow-config",
        str(config_path),
    ]


@pytest.mark.parametrize(
    ("sim_types", "expected"),
    [
        # both PSD methods present (simulate_psd + simulate_psd_with_psl)
        (("single_temp", "pulse_lib"), {"single_template", "psl"}),
        # only the single-template method present (the shipped default, with
        # simulate_psd_with_psl off): the pulse_lib field is absent on disk
        (("single_temp",), {"single_template"}),
    ],
)
def test_extract_hpge_aoemean_energy_dependence_mc_only(
    tmp_path, monkeypatch, sim_types, expected
):
    hit_file = tmp_path / "test-Pb212-hit.lh5"
    _write_fake_mc_hit(hit_file, sim_types=sim_types)

    monkeypatch.setattr(sys, "argv", _build_argv(tmp_path, [hit_file]))

    extract_hpge_aoemean_energy_dependence.main()

    pars_file = tmp_path / "outputs" / "V05261A-model.yaml"
    plot_file = tmp_path / "plots" / "V05261A-fit-result.pdf"
    assert pars_file.exists()
    assert plot_file.exists()

    out = yaml.safe_load(pars_file.read_text())
    assert set(out) == {"energy_corrections", "mc_resolution"}
    assert set(out["energy_corrections"]) == expected
    assert set(out["mc_resolution"]) == expected
    assert "data_validation" not in out

    # the mean energy correction is the primary product: its fit parameters must
    # be real numbers (a NaN would mean the underlying fit silently failed)
    for key in expected:
        pars = out["energy_corrections"][key]["pars"]
        assert np.isfinite(pars["a"])
        assert np.isfinite(pars["b"])
