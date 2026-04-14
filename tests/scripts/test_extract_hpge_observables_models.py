from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from legendsimflow.scripts.pars import extract_hpge_observables_models

dummyprod = Path(__file__).parent.parent / "dummyprod"

# p03 runid: eresmod validity.yaml maps 20230312 → p03 file (has "default" + V02160A override)
RUNID_P03 = "l200-p03-r000-phy"
# p02 runid: eresmod validity.yaml maps 20220102 → p02 file (no "default", only V02160A override)
RUNID_P02 = "l200-p02-r000-phy"

# count of geds detectors in the p03 channelmap (grep "system: geds" gives 101 lines)
N_GEDS_P03 = 101

# fake l200data base returned by the mocked lookup_energy_res_metadata
_FAKE_L200DATA_ERESMOD = {
    "V05261A": {
        "expression": "FWHMLinear",
        "parameters": {"a": 0.7, "b": 0.002},
        "uncertainties": {},
    },
    "V02160A": {
        "expression": "FWHMLinear",
        "parameters": {"a": 0.6, "b": 0.0015},
        "uncertainties": {},
    },
}

# keep RUNID as alias for backward compat within this module
RUNID = RUNID_P03


def _build_argv(
    tmp_path: Path, runid: str = RUNID_P03, l200data: str | None = None
) -> list[str]:
    """Return sys.argv for the script under test."""
    config_path = tmp_path / "simflow-config.yaml"
    raw = yaml.safe_load((dummyprod / "simflow-config.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    if l200data is not None:
        raw["paths"]["l200data"] = l200data
    config_path.write_text(yaml.safe_dump(raw))

    return [
        "extract-hpge-obs-models",
        "--runid",
        runid,
        "--eresmod-file",
        str(tmp_path / "eresmod.yaml"),
        "--aoeresmod-file",
        str(tmp_path / "aoeresmod.yaml"),
        "--psdcuts-file",
        str(tmp_path / "psdcuts.yaml"),
        "--simflow-config",
        str(config_path),
    ]


def _run_with_mocks(
    tmp_path: Path,
    monkeypatch,
    runid: str = RUNID_P03,
    l200data: str | None = None,
) -> Path:
    """Patch sys.argv and the l200data-dependent helpers, then call main().

    The eresmod path is exercised via metadata; the aoeresmod and psdcuts
    paths are short-circuited by returning empty dicts so no l200data is needed.
    Pass a non-None `l200data` to exercise the l200data extraction path (the
    actual extraction is still mocked via lookup_energy_res_metadata).
    """
    monkeypatch.setattr(
        sys, "argv", _build_argv(tmp_path, runid=runid, l200data=l200data)
    )

    with (
        patch(
            "legendsimflow.utils.get_hit_tier_name",
            return_value="hit",
        ),
        patch(
            "legendsimflow.utils.init_generated_pars_db",
            return_value=None,
        ),
        patch(
            "legendsimflow.hpge_pars.lookup_energy_res_metadata",
            return_value=_FAKE_L200DATA_ERESMOD,
        ),
        patch(
            "legendsimflow.hpge_pars.lookup_aoe_res_metadata",
            return_value={},
        ),
        patch(
            "legendsimflow.hpge_pars.lookup_psd_cut_values",
            return_value={},
        ),
    ):
        extract_hpge_observables_models.main()

    return tmp_path / "eresmod.yaml"


def test_extract_eresmod_from_metadata(tmp_path, monkeypatch):
    """Metadata-driven path: eresmod.yaml must be expanded to per-detector entries."""
    eresmod_file = _run_with_mocks(tmp_path, monkeypatch)

    assert eresmod_file.exists(), "eresmod output file was not created"

    with eresmod_file.open() as f:
        result = yaml.safe_load(f)

    assert isinstance(result, dict), "eresmod output must be a YAML mapping"

    # "default" key must NOT appear — non-ON fallback comes from hit tier settings
    assert "default" not in result, "eresmod output must not contain a 'default' key"

    # every per-detector entry must carry expression and parameters with a and b
    for det, entry in result.items():
        assert "expression" in entry, f"{det}: missing 'expression'"
        assert "parameters" in entry, f"{det}: missing 'parameters'"
        assert "a" in entry["parameters"], f"{det}: missing parameter 'a'"
        assert "b" in entry["parameters"], f"{det}: missing parameter 'b'"

    # the per-detector override for V02160A must win over the default
    assert "V02160A" in result, "V02160A must be present in the eresmod output"
    assert result["V02160A"]["parameters"]["a"] == pytest.approx(0.3)
    assert result["V02160A"]["parameters"]["b"] == pytest.approx(0.0009)

    # at least one other geds detector must carry the default values
    default_detectors = [
        det
        for det, entry in result.items()
        if det != "V02160A"
        and entry["parameters"]["a"] == pytest.approx(0.5)
        and entry["parameters"]["b"] == pytest.approx(0.001)
    ]
    assert len(default_detectors) >= 1, (
        "at least one geds detector must have the default a=0.5, b=0.001"
    )


def test_extract_eresmod_per_detector_count(tmp_path, monkeypatch):
    """Metadata-driven path: one output entry per geds detector in the p03 channelmap."""
    eresmod_file = _run_with_mocks(tmp_path, monkeypatch)

    with eresmod_file.open() as f:
        result = yaml.safe_load(f)

    assert len(result) == N_GEDS_P03, (
        f"expected {N_GEDS_P03} geds detectors in eresmod output, got {len(result)}"
    )


def test_extract_eresmod_l200data_with_overrides(tmp_path, monkeypatch):
    """Case 2: l200data base + per-detector overrides from metadata (no 'default' key).

    The p02 eresmod file has only a V02160A override. The l200data base is
    mocked to return V05261A and V02160A. The output must contain both detectors,
    with V02160A replaced by the metadata value and V05261A kept from l200data.
    """
    eresmod_file = _run_with_mocks(
        tmp_path, monkeypatch, runid=RUNID_P02, l200data="/fake/l200data"
    )

    with eresmod_file.open() as f:
        result = yaml.safe_load(f)

    assert isinstance(result, dict), "eresmod output must be a YAML mapping"
    assert "default" not in result, "eresmod output must not contain a 'default' key"

    # V05261A: not in metadata overrides — must keep l200data values
    assert "V05261A" in result, "V05261A must be present (from l200data base)"
    assert result["V05261A"]["parameters"]["a"] == pytest.approx(0.7)
    assert result["V05261A"]["parameters"]["b"] == pytest.approx(0.002)

    # V02160A: in metadata overrides — must be replaced by metadata values
    assert "V02160A" in result, "V02160A must be present (overridden by metadata)"
    assert result["V02160A"]["parameters"]["a"] == pytest.approx(0.3)
    assert result["V02160A"]["parameters"]["b"] == pytest.approx(0.0009)
