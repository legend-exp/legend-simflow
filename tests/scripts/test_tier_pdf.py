from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml
from lgdo import Array, Table, VectorOfVectors, lh5

from legendsimflow.scripts.tier import pdf

dummyprod = Path(__file__).parent.parent / "dummyprod"


def _make_cvt_file(path: Path) -> None:
    """Write a minimal cvt LH5 file with the fields expected by the pdf script.

    6 events:
    - event 0: m1, 500 keV,      spms=False, PSD valid+simulated+single-site → passes all cuts
    - event 1: m1, 1000 keV,     spms=False, PSD valid+simulated+multi-site  → fail/psd
    - event 2: m1, 1500 keV,     spms=False, PSD not valid                   → cut (no observable)
    - event 3: m1, 2000 keV,     spms=True,  PSD valid+simulated+single-site → fail/lar
    - event 4: m2, 200+300 keV,  spms=False
    - event 5: m1, 2500 keV,     spms=False, PSD valid+not simulated         → cut (no observable)
    """
    multiplicity = Array(np.array([1, 1, 1, 1, 2, 1], dtype=np.int32))
    energy = VectorOfVectors(
        data=[[500.0], [1000.0], [1500.0], [2000.0], [200.0, 300.0], [2500.0]]
    )
    is_good_channel = VectorOfVectors(
        data=[[True], [True], [True], [True], [True, True], [True]]
    )
    is_good = VectorOfVectors(
        data=[[True], [True], [False], [True], [True, True], [True]]
    )
    has_aoe = VectorOfVectors(
        data=[[True], [True], [False], [True], [True, True], [False]]
    )
    is_single_site = VectorOfVectors(
        data=[[True], [False], [False], [True], [True, True], [False]]
    )
    spms = Array(np.array([False, False, False, True, False, False]))

    psd = Table(
        col_dict={
            "is_good": is_good,
            "has_aoe": has_aoe,
            "is_single_site": is_single_site,
        }
    )
    geds = Table(
        col_dict={
            "energy": energy,
            "is_good_channel": is_good_channel,
            "multiplicity": multiplicity,
            "psd": psd,
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
    root_keys = lh5.ls(pdf_file)
    assert "pdf" in root_keys
    assert "nr_sim_events" in root_keys
    nr_sim_events = lh5.read("nr_sim_events", pdf_file)
    assert np.issubdtype(type(nr_sim_events.value), np.integer)

    inner_keys = lh5.ls(pdf_file, "pdf/")
    for name in (
        "pdf/hit",
        "pdf/mul",
        "pdf/mul_lar",
        "pdf/mul_psd",
        "pdf/mul_lar_psd",
        "pdf/mul2",
        "pdf/fail",
    ):
        assert name in inner_keys
    fail_keys = lh5.ls(pdf_file, "pdf/fail/")
    for name in ("pdf/fail/lar", "pdf/fail/psd"):
        assert name in fail_keys

    def _sum(path):
        return lh5.read_as(path, pdf_file, "hist").sum()

    # 7 per-detector deposits: 5 m1 + 2 m2
    assert _sum("pdf/hit") == 7

    # 5 m1 events
    assert _sum("pdf/mul") == 5

    # event 3 (spms=True) is vetoed: 4 survive
    assert _sum("pdf/mul_lar") == 4

    # events 0 and 3 pass PSD (valid+simulated+SS); events 1 (MS), 2 (invalid), 5 (not simulated) cut
    assert _sum("pdf/mul_psd") == 2

    # event 3 also fails LAr, so only event 0 survives both
    assert _sum("pdf/mul_lar_psd") == 1

    # 1 m2 event
    assert _sum("pdf/mul2") == 1

    # event 3 fails LAr
    assert _sum("pdf/fail/lar") == 1

    # event 1 (valid+simulated PSD, multi-site) fails PSD; event 5 (not simulated) excluded
    assert _sum("pdf/fail/psd") == 1
