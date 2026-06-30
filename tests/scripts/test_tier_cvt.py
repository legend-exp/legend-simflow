from __future__ import annotations

import sys
from pathlib import Path

import lh5
import numpy as np
import pytest
import yaml
from lgdo import Array, Scalar, Struct, Table

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


def _write_config(tmp_path: Path) -> Path:
    config_path = tmp_path / "simflow-config.yaml"
    raw_config = yaml.safe_load((dummyprod / "simflow-config.yaml").read_text())
    raw_config["paths"]["metadata"] = str(dummyprod / "inputs")
    config_path.write_text(yaml.safe_dump(raw_config))
    return config_path


def _run_cvt(monkeypatch, evt_files: list[Path], cvt_file: Path, config_path: Path):
    argv = [
        "cvt",
        "--evt-files",
        *[str(f) for f in evt_files],
        "--cvt-file",
        str(cvt_file),
        "--simflow-config",
        str(config_path),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    cvt.main()


def test_cvt_script_cli(tmp_path, monkeypatch):
    config_path = _write_config(tmp_path)

    evt_file = tmp_path / "evt.lh5"
    lh5.write(Table(), "evt", evt_file, wo_mode="write_safe")
    lh5.write(Scalar(1000), "number_of_simulated_events", evt_file, wo_mode="append")

    cvt_file = tmp_path / "cvt.lh5"
    _run_cvt(monkeypatch, [evt_file], cvt_file, config_path)

    assert cvt_file.exists()
    assert "evt" in lh5.ls(cvt_file)
    # the single-file copy path must preserve number_of_simulated_events
    assert lh5.read("number_of_simulated_events", cvt_file).value == 1000


def test_cvt_script_cli_sums_number_of_events(tmp_path, monkeypatch):
    """With multiple evt files, number_of_events is summed across jobs."""
    config_path = _write_config(tmp_path)

    evt_files = []
    for i, n in enumerate((1000, 2500)):
        evt_file = tmp_path / f"evt_{i}.lh5"
        lh5.write(
            Table(col_dict={"x": Array(np.arange(3, dtype="int32"))}),
            "evt",
            evt_file,
            wo_mode="write_safe",
        )
        _write_detector_uids(evt_file, {"V01": 11})
        lh5.write(Scalar(n), "number_of_simulated_events", evt_file, wo_mode="append")
        evt_files.append(evt_file)

    cvt_file = tmp_path / "cvt.lh5"
    _run_cvt(monkeypatch, evt_files, cvt_file, config_path)

    assert cvt_file.exists()
    assert "evt" in lh5.ls(cvt_file)
    assert lh5.read("number_of_simulated_events", cvt_file).value == 3500
