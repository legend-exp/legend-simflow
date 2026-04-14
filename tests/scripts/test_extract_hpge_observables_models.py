from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml
from dbetto import AttrsDict

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


_FAKE_L200DATA_AOERESMOD = {
    "V05261A": {
        "expression": "SigmaFit",
        "parameters": {"a": 0.0003, "b": 0, "c": 1},
        "uncertainties": {},
    },
    "V02160A": {
        "expression": "SigmaFit",
        "parameters": {"a": 0.0004, "b": 0, "c": 1},
        "uncertainties": {},
    },
}

_FAKE_L200DATA_PSDCUTS = {
    "V05261A": {"aoe": {"low_side": -2.0, "high_side": 4.0}},
    "V02160A": {"aoe": {"low_side": -1.9, "high_side": 3.9}},
}


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
    fake_aoeresmod: dict | None = None,
    fake_psdcuts: dict | None = None,
) -> tuple[Path, Path, Path]:
    """Patch sys.argv and the l200data-dependent helpers, then call main().

    Returns a tuple of (eresmod_file, aoeresmod_file, psdcuts_file) paths.

    The eresmod path is exercised via metadata; by default aoeresmod and psdcuts
    paths are short-circuited by returning empty dicts so no l200data is needed.
    Pass ``fake_aoeresmod`` or ``fake_psdcuts`` to control what the l200data
    lookup mocks return (used when testing the l200data extraction path).
    Pass a non-None ``l200data`` to exercise the l200data extraction path (the
    actual extraction is still mocked via the lookup_* helpers).
    """
    if fake_aoeresmod is None:
        fake_aoeresmod = {}
    if fake_psdcuts is None:
        fake_psdcuts = {}

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
            return_value=fake_aoeresmod,
        ),
        patch(
            "legendsimflow.hpge_pars.lookup_psd_cut_values",
            return_value=fake_psdcuts,
        ),
    ):
        extract_hpge_observables_models.main()

    return (
        tmp_path / "eresmod.yaml",
        tmp_path / "aoeresmod.yaml",
        tmp_path / "psdcuts.yaml",
    )


def test_extract_eresmod_from_metadata(tmp_path, monkeypatch):
    """Metadata-driven path: eresmod.yaml must be expanded to per-detector entries."""
    eresmod_file, _aoe, _psd = _run_with_mocks(tmp_path, monkeypatch)

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
    eresmod_file, _aoe, _psd = _run_with_mocks(tmp_path, monkeypatch)

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
    eresmod_file, _aoe, _psd = _run_with_mocks(
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


# ---------------------------------------------------------------------------
# A/E resolution model (aoeresmod) tests
# ---------------------------------------------------------------------------


def test_extract_aoeresmod_from_metadata(tmp_path, monkeypatch):
    """Metadata-driven path: aoeresmod.yaml must be expanded to per-detector entries."""
    _eres, aoeresmod_file, _psd = _run_with_mocks(tmp_path, monkeypatch)

    assert aoeresmod_file.exists(), "aoeresmod output file was not created"

    with aoeresmod_file.open() as f:
        result = yaml.safe_load(f)

    assert isinstance(result, dict), "aoeresmod output must be a YAML mapping"

    # "default" key must NOT appear in the expanded per-detector output
    assert "default" not in result, "aoeresmod output must not contain a 'default' key"

    # every per-detector entry must carry expression and parameters with a, b, c
    for det, entry in result.items():
        assert "expression" in entry, f"{det}: missing 'expression'"
        assert "parameters" in entry, f"{det}: missing 'parameters'"
        assert "a" in entry["parameters"], f"{det}: missing parameter 'a'"

    # the per-detector override for V02160A must win over the default (a=0.0001)
    assert "V02160A" in result, "V02160A must be present in the aoeresmod output"
    assert result["V02160A"]["parameters"]["a"] == pytest.approx(0.0002)

    # at least one other geds detector must carry the default values (a=0.0001)
    default_detectors = [
        det
        for det, entry in result.items()
        if det != "V02160A" and entry["parameters"]["a"] == pytest.approx(0.0001)
    ]
    assert len(default_detectors) >= 1, (
        "at least one geds detector must have the default a=0.0001"
    )


def test_extract_aoeresmod_per_detector_count(tmp_path, monkeypatch):
    """Metadata-driven path: one output entry per geds detector in the p03 channelmap."""
    _eres, aoeresmod_file, _psd = _run_with_mocks(tmp_path, monkeypatch)

    with aoeresmod_file.open() as f:
        result = yaml.safe_load(f)

    assert len(result) == N_GEDS_P03, (
        f"expected {N_GEDS_P03} geds detectors in aoeresmod output, got {len(result)}"
    )


def test_extract_aoeresmod_l200data_with_overrides(tmp_path, monkeypatch):
    """Case 2: l200data base + per-detector overrides from metadata (no 'default' key).

    The p02 aoeresmod file has only a V02160A override. The l200data base is
    mocked to return V05261A and V02160A. The output must contain both detectors,
    with V02160A replaced by the metadata value (a=0.0002) and V05261A kept
    from l200data (a=0.0003).
    """
    _, aoeresmod_file, _ = _run_with_mocks(
        tmp_path,
        monkeypatch,
        runid=RUNID_P02,
        l200data="/fake/l200data",
        fake_aoeresmod=_FAKE_L200DATA_AOERESMOD,
    )

    with aoeresmod_file.open() as f:
        result = yaml.safe_load(f)

    assert isinstance(result, dict), "aoeresmod output must be a YAML mapping"
    assert "default" not in result, "aoeresmod output must not contain a 'default' key"

    # V05261A: not in metadata overrides — must keep l200data values (a=0.0003)
    assert "V05261A" in result, "V05261A must be present (from l200data base)"
    assert result["V05261A"]["parameters"]["a"] == pytest.approx(0.0003)

    # V02160A: in metadata overrides — must be replaced by metadata values (a=0.0002)
    assert "V02160A" in result, "V02160A must be present (overridden by metadata)"
    assert result["V02160A"]["parameters"]["a"] == pytest.approx(0.0002)


# ---------------------------------------------------------------------------
# PSD cut values (psdcuts) tests
# ---------------------------------------------------------------------------


def test_extract_psdcuts_from_metadata(tmp_path, monkeypatch):
    """Metadata-driven path: psdcuts.yaml must be expanded to per-detector entries."""
    _eres, _aoe, psdcuts_file = _run_with_mocks(tmp_path, monkeypatch)

    assert psdcuts_file.exists(), "psdcuts output file was not created"

    with psdcuts_file.open() as f:
        result = yaml.safe_load(f)

    assert isinstance(result, dict), "psdcuts output must be a YAML mapping"

    # "default" key must NOT appear in the expanded per-detector output
    assert "default" not in result, "psdcuts output must not contain a 'default' key"

    # every per-detector entry must carry an aoe sub-dict with low_side and high_side
    for det, entry in result.items():
        assert "aoe" in entry, f"{det}: missing 'aoe' key"
        assert "low_side" in entry["aoe"], f"{det}: missing 'aoe.low_side'"
        assert "high_side" in entry["aoe"], f"{det}: missing 'aoe.high_side'"

    # the per-detector override for V02160A must win over the default (low_side=-1.5)
    assert "V02160A" in result, "V02160A must be present in the psdcuts output"
    assert result["V02160A"]["aoe"]["low_side"] == pytest.approx(-1.4)

    # at least one other geds detector must carry the default values (low_side=-1.5)
    default_detectors = [
        det
        for det, entry in result.items()
        if det != "V02160A" and entry["aoe"]["low_side"] == pytest.approx(-1.5)
    ]
    assert len(default_detectors) >= 1, (
        "at least one geds detector must have the default aoe.low_side=-1.5"
    )


def test_extract_psdcuts_per_detector_count(tmp_path, monkeypatch):
    """Metadata-driven path: one output entry per geds detector in the p03 channelmap."""
    _eres, _aoe, psdcuts_file = _run_with_mocks(tmp_path, monkeypatch)

    with psdcuts_file.open() as f:
        result = yaml.safe_load(f)

    assert len(result) == N_GEDS_P03, (
        f"expected {N_GEDS_P03} geds detectors in psdcuts output, got {len(result)}"
    )


def test_extract_psdcuts_l200data_with_overrides(tmp_path, monkeypatch):
    """Case 2: l200data base + per-detector overrides from metadata (no 'default' key).

    The p02 psdcuts file has only a V02160A override. The l200data base is
    mocked to return V05261A and V02160A. The output must contain both detectors,
    with V02160A replaced by the metadata value (low_side=-1.4) and V05261A kept
    from l200data (low_side=-2.0).
    """
    _, _, psdcuts_file = _run_with_mocks(
        tmp_path,
        monkeypatch,
        runid=RUNID_P02,
        l200data="/fake/l200data",
        fake_psdcuts=_FAKE_L200DATA_PSDCUTS,
    )

    with psdcuts_file.open() as f:
        result = yaml.safe_load(f)

    assert isinstance(result, dict), "psdcuts output must be a YAML mapping"
    assert "default" not in result, "psdcuts output must not contain a 'default' key"

    # V05261A: not in metadata overrides — must keep l200data values (low_side=-2.0)
    assert "V05261A" in result, "V05261A must be present (from l200data base)"
    assert result["V05261A"]["aoe"]["low_side"] == pytest.approx(-2.0)

    # V02160A: in metadata overrides — must be replaced by metadata values (low_side=-1.4)
    assert "V02160A" in result, "V02160A must be present (overridden by metadata)"
    assert result["V02160A"]["aoe"]["low_side"] == pytest.approx(-1.4)


# ---------------------------------------------------------------------------
# Error-path tests: no l200data and no metadata default
# ---------------------------------------------------------------------------


def _make_eresmod_default():
    return AttrsDict(
        {
            "default": AttrsDict(
                {
                    "expression": "FWHMLinear",
                    "parameters": AttrsDict({"a": 0.5, "b": 0.001}),
                }
            )
        }
    )


def _make_aoeresmod_default():
    return AttrsDict(
        {
            "default": AttrsDict(
                {
                    "expression": "SigmaFit",
                    "parameters": AttrsDict({"a": 0.0001, "b": 0, "c": 1}),
                }
            )
        }
    )


def test_extract_aoeresmod_raises_without_l200data_and_default(tmp_path, monkeypatch):
    """RuntimeError when l200data is absent and no aoeresmod 'default' key exists."""
    monkeypatch.setattr(
        sys, "argv", _build_argv(tmp_path, runid=RUNID_P03, l200data=None)
    )

    def _simpars(_metadata, par, _runid, **_kw):
        if par == "geds.eresmod":
            return _make_eresmod_default()
        return None

    with (
        patch("legendsimflow.metadata.simpars", side_effect=_simpars),
        pytest.raises(RuntimeError, match="aoeresmod"),
    ):
        extract_hpge_observables_models.main()


def test_extract_psdcuts_raises_without_l200data_and_default(tmp_path, monkeypatch):
    """RuntimeError when l200data is absent and no psdcuts 'default' key exists."""
    monkeypatch.setattr(
        sys, "argv", _build_argv(tmp_path, runid=RUNID_P03, l200data=None)
    )

    def _simpars(_metadata, par, _runid, **_kw):
        if par == "geds.eresmod":
            return _make_eresmod_default()
        if par == "geds.aoeresmod":
            return _make_aoeresmod_default()
        return None

    with (
        patch("legendsimflow.metadata.simpars", side_effect=_simpars),
        pytest.raises(RuntimeError, match="psdcuts"),
    ):
        extract_hpge_observables_models.main()
