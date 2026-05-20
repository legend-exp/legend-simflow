from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

from legendsimflow.scripts import extract_hpge_current_pulse_model

dummyprod = Path(__file__).parent.parent / "dummyprod"
l200data = Path(__file__).parent.parent / "l200data" / "v3.0.0"

# p03 runid: currmod validity.yaml maps 20230312 → p03 file (has "default" + V02160A override)
RUNID_P03 = "l200-p03-r000-phy"
# p02 runid: currmod validity.yaml maps 20220102 → p02 file (no "default", only V02160A override)
RUNID_P02 = "l200-p02-r000-phy"

# any ON detector that is NOT V02160A — receives the metadata default values
DET_DEFAULT = "V05261A"
# has a per-detector override in the p03 metadata
DET_OVERRIDE = "V02160A"


def _build_argv(
    tmp_path: Path,
    runid: str,
    hpge_detector: str,
    l200data_path: str | None = None,
) -> list[str]:
    """Return sys.argv for the currmod script under test.

    Copies the dummyprod simflow-config.yaml to tmp_path, points metadata at
    the dummyprod inputs, and ensures no l200data key is present.
    """
    config_path = tmp_path / "simflow-config.yaml"
    raw = yaml.safe_load((dummyprod / "simflow-config.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    if l200data_path is None:
        raw["paths"].pop("l200data", None)
    else:
        raw["paths"]["l200data"] = l200data_path
    config_path.write_text(yaml.safe_dump(raw))

    return [
        "extract-hpge-currmod",
        "--runid",
        runid,
        "--hpge-detector",
        hpge_detector,
        "--pars-file",
        str(tmp_path / "pars.yaml"),
        "--plot-file",
        str(tmp_path / "plot.pdf"),
        "--simflow-config",
        str(config_path),
    ]


def test_metadata_default_written_for_detector(tmp_path, monkeypatch):
    """p03 + a detector without its own override must get the metadata default."""
    decorate_calls = []

    monkeypatch.setattr(
        extract_hpge_current_pulse_model,
        "decorate",
        lambda fig: decorate_calls.append(fig),
    )
    monkeypatch.setattr(sys, "argv", _build_argv(tmp_path, RUNID_P03, DET_DEFAULT))

    extract_hpge_current_pulse_model.main()

    pars_file = tmp_path / "pars.yaml"
    plot_file = tmp_path / "plot.pdf"

    assert pars_file.exists(), "pars output file was not created"
    assert plot_file.exists(), "plot output file was not created"

    with pars_file.open() as f:
        result = yaml.safe_load(f)

    assert "current_pulse_pars" in result
    assert result["current_pulse_pars"]["sigma"] == pytest.approx(45.0)
    assert result["mean_aoe"] == pytest.approx(0.72)
    assert result["current_reso"] == pytest.approx(4.1)
    assert len(decorate_calls) == 1


def test_metadata_override_wins_for_v02160a(tmp_path, monkeypatch):
    """p03 + V02160A must get the per-detector override, not the default."""
    monkeypatch.setattr(sys, "argv", _build_argv(tmp_path, RUNID_P03, DET_OVERRIDE))

    extract_hpge_current_pulse_model.main()

    pars_file = tmp_path / "pars.yaml"
    assert pars_file.exists(), "pars output file was not created"

    with pars_file.open() as f:
        result = yaml.safe_load(f)

    assert result["current_pulse_pars"]["sigma"] == pytest.approx(48.0)
    assert result["current_reso"] == pytest.approx(5.3)


def test_raises_without_default_and_no_l200data(tmp_path, monkeypatch):
    """p02 has no 'default' key and no l200data is configured → RuntimeError."""
    monkeypatch.setattr(sys, "argv", _build_argv(tmp_path, RUNID_P02, DET_DEFAULT))

    with pytest.raises(RuntimeError, match=r"l200data|currmod"):
        extract_hpge_current_pulse_model.main()


def test_l200data_path_creates_valid_outputs(tmp_path, monkeypatch, test_make_ssc_data):
    monkeypatch.setattr(
        sys,
        "argv",
        _build_argv(
            tmp_path,
            runid="l200-p16-r008-ssc",
            hpge_detector="V03422A",
            l200data_path=str(l200data),
        ),
    )

    extract_hpge_current_pulse_model.main()

    pars_file = tmp_path / "pars.yaml"
    plot_file = tmp_path / "plot.pdf"

    assert pars_file.exists(), "pars output file was not created"
    assert plot_file.exists(), "plot output file was not created"

    with pars_file.open() as f:
        result = yaml.safe_load(f)

    assert isinstance(result, dict), "pars output must be a YAML mapping"
    assert "current_pulse_pars" in result
    assert "mean_aoe" in result
    assert "current_reso" in result

    for field in ("height", "t0", "s1", "r1", "w2", "r2", "w3"):
        assert field in result["current_pulse_pars"], (
            f"missing current_pulse_pars.{field}"
        )

    assert isinstance(result["mean_aoe"], (int, float))
    assert isinstance(result["current_reso"], (int, float))
