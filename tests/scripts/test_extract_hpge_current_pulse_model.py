from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

from legendsimflow.scripts import extract_hpge_current_pulse_model

dummyprod = Path(__file__).parent.parent / "dummyprod"

# p03 runid: currmod validity.yaml maps 20230312 → p03 file (has "default" + V02160A override)
RUNID_P03 = "l200-p03-r000-phy"
# p02 runid: currmod validity.yaml maps 20220102 → p02 file (no "default", only V02160A override)
RUNID_P02 = "l200-p02-r000-phy"

# ON detector that receives the metadata default values
DET_DEFAULT = "V05261B"
# has a per-detector override in the p03 metadata
DET_OVERRIDE = "V02160A"


def _build_argv(tmp_path: Path, runid: str) -> list[str]:
    """Return sys.argv for the currmod script under test.

    Copies the dummyprod simflow-config.yaml to tmp_path, points metadata at
    the dummyprod inputs, and ensures no l200data key is present.
    """
    config_path = tmp_path / "simflow-config.yaml"
    raw = yaml.safe_load((dummyprod / "simflow-config.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    raw["paths"].pop("l200data", None)
    config_path.write_text(yaml.safe_dump(raw))

    # minimal cache file with the two dummyprod ON detectors
    cache_path = tmp_path / "modelable_hpge_detectors.yaml"
    cache_path.write_text(
        yaml.safe_dump({runid: {DET_DEFAULT: 4200, DET_OVERRIDE: 4200}})
    )

    return [
        "extract-hpge-currmod",
        "--runid",
        runid,
        "--modelable-hpges-file",
        str(cache_path),
        "--pars-file",
        str(tmp_path / "pars.yaml"),
        "--plot-file",
        str(tmp_path / "plot.pdf"),
        "--simflow-config",
        str(config_path),
    ]


def test_metadata_all_detectors(tmp_path, monkeypatch):
    """p03 with metadata defaults: output is keyed by detector name."""
    monkeypatch.setattr(sys, "argv", _build_argv(tmp_path, RUNID_P03))

    extract_hpge_current_pulse_model.main()

    pars_file = tmp_path / "pars.yaml"
    plot_file = tmp_path / "plot.pdf"

    assert pars_file.exists(), "pars output file was not created"
    assert plot_file.exists(), "plot output file was not created"

    with pars_file.open() as f:
        result = yaml.safe_load(f)

    # DET_DEFAULT (V05261B) gets the metadata default values
    assert DET_DEFAULT in result
    assert result[DET_DEFAULT]["current_pulse_pars"]["sigma"] == pytest.approx(0.1)
    assert result[DET_DEFAULT]["mean_aoe"] == pytest.approx(1.0)
    assert result[DET_DEFAULT]["current_reso"] == pytest.approx(0.01)

    # DET_OVERRIDE (V02160A) gets the per-detector override
    assert DET_OVERRIDE in result
    assert result[DET_OVERRIDE]["current_pulse_pars"]["sigma"] == pytest.approx(0.12)
    assert result[DET_OVERRIDE]["current_reso"] == pytest.approx(0.012)


def test_raises_without_default_and_no_l200data(tmp_path, monkeypatch):
    """p02 has no 'default' key and no l200data is configured → RuntimeError."""
    monkeypatch.setattr(sys, "argv", _build_argv(tmp_path, RUNID_P02))

    with pytest.raises(RuntimeError, match=r"l200data|currmod"):
        extract_hpge_current_pulse_model.main()
