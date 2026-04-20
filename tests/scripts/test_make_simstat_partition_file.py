from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

from legendsimflow.scripts.make_simstat_partition_file import main

dummyprod = Path(__file__).parent.parent / "dummyprod"


@pytest.mark.needs_remage
def test_simstat_partition_file(tmp_path, legend_stp_path, monkeypatch):
    config_path = tmp_path / "simflow-config-l1000.yaml"
    raw = yaml.safe_load((dummyprod / "simflow-config-l1000.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    raw["paths"]["generated"] = str(tmp_path)
    raw["paths"]["pars"] = str(tmp_path / "pars")
    config_path.write_text(yaml.safe_dump(raw))

    output_file = tmp_path / "partitions_ultem_insulators_Pb214_to_Po214.yaml"

    argv = [
        "make-simstat-partition-file",
        "--stp-files",
        str(legend_stp_path),
        "--runlist",
        "l200-p03-r000-phy",
        "l200-p03-r001-phy",
        "--output-file",
        str(output_file),
        "--simflow-config",
        str(config_path),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    main()

    assert output_file.exists()

    result = yaml.safe_load(output_file.read_text())
    assert "job_0000" in result

    job = result["job_0000"]
    assert len(job) > 0

    # all event ranges must be non-negative and ordered
    for _, (start, end) in job.items():
        assert start >= 0
        assert end >= start
