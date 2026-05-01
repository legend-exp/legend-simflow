from __future__ import annotations

import sys
from pathlib import Path

import lh5
import yaml
from lgdo import Table

from legendsimflow.scripts.tier import cvt

dummyprod = Path(__file__).parent.parent / "dummyprod"


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
