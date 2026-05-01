from __future__ import annotations

import sys
from pathlib import Path

import lh5
import pytest
import yaml
from lgdo import Scalar, Struct, Table

from legendsimflow.scripts.tier import cvt

dummyprod = Path(__file__).parent.parent / "dummyprod"


def _write_detector_uids(path: Path, mapping: dict[str, int]) -> None:
    s = Struct({name: Scalar(uid) for name, uid in mapping.items()})
    lh5.write(s, "detector_uids", path, wo_mode="write_safe")


def test_union_detector_uids_identical_mappings(tmp_path):
    f1 = tmp_path / "a.lh5"
    f2 = tmp_path / "b.lh5"
    _write_detector_uids(f1, {"V01": 11, "V02": 12})
    _write_detector_uids(f2, {"V01": 11, "V02": 12})

    assert cvt.union_detector_uids([f1, f2]) == {"V01": 11, "V02": 12}


def test_union_detector_uids_disjoint_subsets(tmp_path):
    """A detector present in one job but missing in another (zero-hit) must union."""
    f1 = tmp_path / "a.lh5"
    f2 = tmp_path / "b.lh5"
    _write_detector_uids(f1, {"V01": 11, "V02": 12})
    _write_detector_uids(f2, {"V01": 11, "V03": 13})

    assert cvt.union_detector_uids([f1, f2]) == {"V01": 11, "V02": 12, "V03": 13}


def test_union_detector_uids_name_collision_raises(tmp_path):
    f1 = tmp_path / "a.lh5"
    f2 = tmp_path / "b.lh5"
    _write_detector_uids(f1, {"V01": 11})
    _write_detector_uids(f2, {"V01": 99})

    with pytest.raises(ValueError, match="'V01' maps to 11"):
        cvt.union_detector_uids([f1, f2])


def test_union_detector_uids_uid_collision_raises(tmp_path):
    f1 = tmp_path / "a.lh5"
    f2 = tmp_path / "b.lh5"
    _write_detector_uids(f1, {"V01": 11})
    _write_detector_uids(f2, {"V02": 11})

    with pytest.raises(ValueError, match="uid 11 maps to"):
        cvt.union_detector_uids([f1, f2])


def test_cvt_script_cli(tmp_path, monkeypatch):
    config_path = tmp_path / "simflow-config.yaml"
    raw_config = yaml.safe_load((dummyprod / "simflow-config.yaml").read_text())
    raw_config["paths"]["metadata"] = str(dummyprod / "inputs")
    config_path.write_text(yaml.safe_dump(raw_config))

    evt_file = tmp_path / "evt.lh5"
    lh5.write(Table(), "evt", evt_file, wo_mode="write_safe")

    cvt_file = tmp_path / "cvt.lh5"

    argv = [
        "cvt",
        "--evt-files",
        str(evt_file),
        "--cvt-file",
        str(cvt_file),
        "--simflow-config",
        str(config_path),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    cvt.main()

    assert cvt_file.exists()
    assert "evt" in lh5.ls(cvt_file)
