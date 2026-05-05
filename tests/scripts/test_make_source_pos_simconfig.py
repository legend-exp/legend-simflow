from __future__ import annotations

import copy
import sys
from pathlib import Path

import dbetto.utils
import pytest
import yaml

from legendsimflow.scripts.make_source_pos_simconfig import (
    JOBS,
    PRIMARIES,
    _swap_tier_segment,
    main,
    replace_position,
)

# ---------------------------------------------------------------------------
# _swap_tier_segment
# ---------------------------------------------------------------------------


def test_swap_tier_segment_basic():
    result = _swap_tier_segment("/some/path/tier/stp/experiment", "stp", "hit")
    assert result == "/some/path/tier/hit/experiment"


def test_swap_tier_segment_only_last_occurrence():
    """A parent directory whose name contains the tier string must not be touched."""
    result = _swap_tier_segment("/stp-study/config/tier/stp/exp", "stp", "hit")
    assert result == "/stp-study/config/tier/hit/exp"


def test_swap_tier_segment_missing_raises():
    with pytest.raises(ValueError, match="not found as a path component"):
        _swap_tier_segment("/some/path/hit/experiment", "stp", "hit")


# ---------------------------------------------------------------------------
# replace_position
# ---------------------------------------------------------------------------

_VALID_CONFIG = {
    "geom_config_extra": {
        "sis": {
            "1": {
                "sis_z": 8224,
                "offset": -63,
                "sources": [None, None, "Th228", None],
            }
        }
    },
    "primaries_per_job": 100_000,
    "number_of_jobs": 1,
}


def test_replace_position_sets_offsets():
    cfg = copy.deepcopy(_VALID_CONFIG)
    result = replace_position(cfg, dz=10, dphi=90)
    sis_cfg = result["geom_config_extra"]["sis"]["1"]
    assert sis_cfg["offset"] == 10
    assert sis_cfg["phi_offset"] == 90


def test_replace_position_sets_primaries_and_jobs():
    cfg = copy.deepcopy(_VALID_CONFIG)
    result = replace_position(cfg, dz=0, dphi=0)
    assert result["primaries_per_job"] == PRIMARIES
    assert result["number_of_jobs"] == JOBS


def test_replace_position_missing_geom_config_extra():
    with pytest.raises(ValueError, match="geom_config_extra"):
        replace_position({}, dz=0, dphi=0)


def test_replace_position_missing_sis():
    with pytest.raises(ValueError, match="sis"):
        replace_position({"geom_config_extra": {}}, dz=0, dphi=0)


def test_replace_position_multiple_sis():
    cfg = copy.deepcopy(_VALID_CONFIG)
    cfg["geom_config_extra"]["sis"]["2"] = {"offset": 0}
    with pytest.raises(ValueError, match="one SIS"):
        replace_position(cfg, dz=0, dphi=0)


# ---------------------------------------------------------------------------
# main() CLI
# ---------------------------------------------------------------------------


def _make_minimal_stp_dir(base: Path, experiment: str) -> Path:
    """Create a minimal stp config directory under base/tier/stp/<experiment>."""
    stp_dir = base / "tier" / "stp" / experiment
    stp_dir.mkdir(parents=True)
    simconfig = {
        "sis_src_Ra224_to_Pb208": {
            "template": "$_/template.mac",
            "geom_config_extra": {
                "sis": {
                    "1": {
                        "sis_z": 8224,
                        "offset": -63,
                        "sources": [None, None, "Th228", None],
                    }
                }
            },
            "generator": "~defines:Ra224_to_Pb208",
            "primaries_per_job": 10_000,
            "number_of_jobs": 2,
        }
    }
    (stp_dir / "simconfig.yaml").write_text(yaml.dump(simconfig))
    return stp_dir


def _make_minimal_hit_dir(base: Path, experiment: str) -> Path:
    """Create a matching hit config directory under base/tier/hit/<experiment>."""
    hit_dir = base / "tier" / "hit" / experiment
    hit_dir.mkdir(parents=True)
    hit_simconfig = {
        "sis_src_Ra224_to_Pb208": {
            "runlist": ["l200-p02-r000-phy"],
        }
    }
    (hit_dir / "simconfig.yaml").write_text(yaml.dump(hit_simconfig))
    return hit_dir


def test_main_creates_expanded_simconfigs(tmp_path, monkeypatch):
    """main() must write expanded stp and hit simconfigs over the z x phi grid."""
    exp = "testexp"
    stp_in = _make_minimal_stp_dir(tmp_path / "input", exp)
    _make_minimal_hit_dir(tmp_path / "input", exp)

    stp_out = tmp_path / "output" / "tier" / "stp" / exp

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "make-source-pos-simconfig",
            "--input-path",
            str(stp_in),
            "--output-path",
            str(stp_out),
            "--z",
            "0",
            "5",
            "--phi",
            "0",
            "90",
        ],
    )
    main()

    # stp simconfig
    stp_simconfig = dbetto.utils.load_dict(stp_out / "simconfig.yaml")
    assert "sis_src_Ra224_to_Pb208_z_0_mm_phi_0_deg" in stp_simconfig
    assert "sis_src_Ra224_to_Pb208_z_0_mm_phi_90_deg" in stp_simconfig
    assert "sis_src_Ra224_to_Pb208_z_5_mm_phi_0_deg" in stp_simconfig
    assert "sis_src_Ra224_to_Pb208_z_5_mm_phi_90_deg" in stp_simconfig
    assert len(stp_simconfig) == 4  # 2 z x 2 phi

    # offsets are set correctly
    entry = stp_simconfig["sis_src_Ra224_to_Pb208_z_5_mm_phi_90_deg"]
    assert entry["geom_config_extra"]["sis"]["1"]["offset"] == 5
    assert entry["geom_config_extra"]["sis"]["1"]["phi_offset"] == 90

    # hit simconfig
    hit_out = tmp_path / "output" / "tier" / "hit" / exp
    hit_simconfig = dbetto.utils.load_dict(hit_out / "simconfig.yaml")
    assert "sis_src_Ra224_to_Pb208_z_0_mm_phi_0_deg" in hit_simconfig
    assert len(hit_simconfig) == 4


def test_main_rejects_same_input_output(tmp_path, monkeypatch):
    exp = "testexp"
    stp_in = _make_minimal_stp_dir(tmp_path / "input", exp)
    _make_minimal_hit_dir(tmp_path / "input", exp)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "make-source-pos-simconfig",
            "--input-path",
            str(stp_in),
            "--output-path",
            str(stp_in),
            "--z",
            "0",
            "--phi",
            "0",
        ],
    )
    with pytest.raises(ValueError, match="same"):
        main()


def test_main_rejects_non_integer_z(tmp_path, monkeypatch):
    exp = "testexp"
    stp_in = _make_minimal_stp_dir(tmp_path / "input", exp)
    _make_minimal_hit_dir(tmp_path / "input", exp)
    stp_out = tmp_path / "output" / "tier" / "stp" / exp

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "make-source-pos-simconfig",
            "--input-path",
            str(stp_in),
            "--output-path",
            str(stp_out),
            "--z",
            "0.5",
            "--phi",
            "0",
        ],
    )
    with pytest.raises(ValueError, match="not an integer"):
        main()


def test_main_rejects_non_integer_phi(tmp_path, monkeypatch):
    exp = "testexp"
    stp_in = _make_minimal_stp_dir(tmp_path / "input", exp)
    _make_minimal_hit_dir(tmp_path / "input", exp)
    stp_out = tmp_path / "output" / "tier" / "stp" / exp

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "make-source-pos-simconfig",
            "--input-path",
            str(stp_in),
            "--output-path",
            str(stp_out),
            "--z",
            "0",
            "--phi",
            "45.5",
        ],
    )
    with pytest.raises(ValueError, match="not an integer"):
        main()
