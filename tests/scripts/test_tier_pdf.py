from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml
from lgdo import Array, Table, VectorOfVectors, lh5

from legendsimflow.scripts.tier import pdf

dummyprod = Path(__file__).parent.parent / "dummyprod"


def _make_cvt_file(path: Path) -> None:
    """Write a minimal cvt LH5 file with the fields expected by the pdf script."""
    # 3 events: multiplicities 1, 1, 2
    multiplicity = Array(np.array([1, 1, 2], dtype=np.int32))
    # per-event energy lists (VectorOfVectors): one or two detectors
    energy = VectorOfVectors(data=[[500.0], [1000.0], [200.0, 300.0]])
    is_good_channel = VectorOfVectors(data=[[True], [True], [True, True]])
    spms = Array(np.array([False, False, False]))

    geds = Table(
        col_dict={
            "energy": energy,
            "is_good_channel": is_good_channel,
            "multiplicity": multiplicity,
        }
    )
    coincident = Table(col_dict={"spms": spms})
    evt = Table(col_dict={"geds": geds, "coincident": coincident})
    lh5.write(evt, "evt", path, wo_mode="write_safe")


def test_pdf_script_cli(tmp_path, monkeypatch):
    config_path = tmp_path / "simflow-config.yaml"
    raw_config = yaml.safe_load((dummyprod / "simflow-config.yaml").read_text())
    raw_config["paths"]["metadata"] = str(dummyprod / "inputs")
    config_path.write_text(yaml.safe_dump(raw_config))

    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)

    pdf_file = tmp_path / "pdf.lh5"

    # use a simid that exists in the dummyprod stp simconfig
    simid = "birds_nest_K40"

    argv = [
        "pdf",
        "--cvt-file",
        str(cvt_file),
        "--pdf-file",
        str(pdf_file),
        "--simid",
        simid,
        "--simflow-config",
        str(config_path),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    pdf.main()

    assert pdf_file.exists()
    assert "pdf" in lh5.ls(pdf_file)
    inner_keys = lh5.ls(pdf_file, "pdf/")
    assert "pdf/n_prim" in inner_keys
    assert "pdf/pdfs" in inner_keys
