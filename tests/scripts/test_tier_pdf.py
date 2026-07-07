from __future__ import annotations

import shutil
import sys
from pathlib import Path

import lh5
import numpy as np
import pytest
import yaml
from lgdo import Array, Scalar, Struct, Table, VectorOfVectors

from legendsimflow.scripts.tier import pdf

dummyprod = Path(__file__).parent.parent / "dummyprod"

# number of simulated primary events stored in the synthetic cvt files; the pdf
# script reads it back and writes it as ``nr_sim_events``.
_N_SIM_EVENTS = 1000


def _write_cvt_root(
    path: Path, names: list[str], uids: list[int], n_sim_events: int = _N_SIM_EVENTS
) -> None:
    """Write the cvt root metadata: detector_uids and number_of_simulated_events."""
    detector_uids = Struct(
        {name: Scalar(int(uid)) for name, uid in zip(names, uids, strict=True)}
    )
    lh5.write(detector_uids, "detector_uids", path, wo_mode="append")
    lh5.write(
        Scalar(n_sim_events), "number_of_simulated_events", path, wo_mode="append"
    )


_SIMID_LEGEND = "birds_nest_K40"
_SIMID_L1000 = "ultem_insulators_Pb214_to_Po214"


def _make_metadata_with_pdf_settings(tmp_path: Path, extra_pdf_settings: dict) -> Path:
    """Copy dummyprod inputs to tmp_path and patch the pdf/legend/settings.yaml."""
    meta_dir = tmp_path / "inputs"
    shutil.copytree(dummyprod / "inputs", meta_dir)
    settings_path = (
        meta_dir / "simprod" / "config" / "tier" / "pdf" / "legend" / "settings.yaml"
    )
    base = yaml.safe_load(settings_path.read_text())
    base.update(extra_pdf_settings)
    settings_path.write_text(yaml.safe_dump(base))
    return meta_dir


def _run_pdf(
    tmp_path: Path,
    monkeypatch,
    cvt_file: Path,
    meta_dir: Path | None = None,
    simid: str = _SIMID_LEGEND,
    config_template: str = "simflow-config.yaml",
) -> Path:
    config_path = tmp_path / config_template
    raw_config = yaml.safe_load((dummyprod / config_template).read_text())
    raw_config["paths"]["metadata"] = str(meta_dir or (dummyprod / "inputs"))
    config_path.write_text(yaml.safe_dump(raw_config))

    pdf_file = tmp_path / "pdf.lh5"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "pdf",
            "--cvt-file",
            str(cvt_file),
            "--pdf-file",
            str(pdf_file),
            "--simid",
            simid,
            "--simflow-config",
            str(config_path),
        ],
    )
    pdf.main()
    return pdf_file


def _make_cvt_file(path: Path) -> None:
    """Write a minimal cvt LH5 file with the fields expected by the pdf script."""
    # 6 events:
    # - event 0: m1, 500 keV,      spms=False, PSD valid+simulated+single-site → passes all cuts
    # - event 1: m1, 1000 keV,     spms=False, PSD valid+simulated+multi-site  → fail/psd
    # - event 2: m1, 1500 keV,     spms=False, PSD not valid                   → cut (no observable)
    # - event 3: m1, 2000 keV,     spms=True,  PSD valid+simulated+single-site → fail/lar
    # - event 4: m2, 200+300 keV,  spms=False
    # - event 5: m1, 2500 keV,     spms=False, PSD valid+not simulated         → cut (no observable)
    multiplicity = Array(np.array([1, 1, 1, 1, 2, 1], dtype=np.int32))
    energy = VectorOfVectors(
        data=[[500.0], [1000.0], [1500.0], [2000.0], [200.0, 300.0], [2500.0]]
    )
    rawid = VectorOfVectors(data=[[1], [1], [1], [1], [1, 2], [1]])
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

    single_temp = Table(
        col_dict={
            "has_aoe": has_aoe,
            "is_single_site": is_single_site,
        }
    )
    psd = Table(
        col_dict={
            "is_good": is_good,
            "single_temp": single_temp,
        }
    )
    geds = Table(
        col_dict={
            "energy": energy,
            "rawid": rawid,
            "is_good_channel": is_good_channel,
            "multiplicity": multiplicity,
            "psd": psd,
        }
    )
    coincident = Table(col_dict={"spms": spms})
    evt = Table(col_dict={"geds": geds, "coincident": coincident})
    lh5.write(evt, "evt", path, wo_mode="write_safe")
    _write_cvt_root(path, ["V01", "B02"], [1, 2])


def test_pdf_script_cli(tmp_path, monkeypatch):
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)
    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file)

    assert pdf_file.exists()
    root_keys = lh5.ls(pdf_file)
    assert "pdf" in root_keys
    assert "nr_sim_events" in root_keys
    nr_sim_events = lh5.read("nr_sim_events", pdf_file)
    assert np.issubdtype(type(nr_sim_events.value), np.integer)
    # the pdf script must read number_of_events back from the cvt file
    assert nr_sim_events.value == _N_SIM_EVENTS

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
    assert _sum("pdf/hit/all") == 7

    # 5 m1 events
    assert _sum("pdf/mul/all") == 5

    # event 3 (spms=True) is vetoed: 4 survive
    assert _sum("pdf/mul_lar/all") == 4

    # events 0 and 3 pass PSD (valid+simulated+SS); events 1 (MS), 2 (invalid), 5 (not simulated) cut
    assert _sum("pdf/mul_psd/all") == 2

    # event 3 also fails LAr, so only event 0 survives both
    assert _sum("pdf/mul_lar_psd/all") == 1

    # 1 m2 event
    assert _sum("pdf/mul2") == 1

    # event 3 fails LAr
    assert _sum("pdf/fail/lar/all") == 1

    # event 1 (valid+simulated PSD, multi-site) fails PSD; event 5 (not simulated) excluded
    assert _sum("pdf/fail/psd/all") == 1


def _make_cvt_file_no_spms(path: Path) -> None:
    """Write a minimal cvt LH5 file with no spms field (mimics ``--skip-opt``)."""
    multiplicity = Array(np.array([1, 1, 1, 1, 2, 1], dtype=np.int32))
    energy = VectorOfVectors(
        data=[[500.0], [1000.0], [1500.0], [2000.0], [200.0, 300.0], [2500.0]]
    )
    rawid = VectorOfVectors(data=[[1], [1], [1], [1], [1, 2], [1]])
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
    single_temp = Table(
        col_dict={
            "has_aoe": has_aoe,
            "is_single_site": is_single_site,
        }
    )
    psd = Table(
        col_dict={
            "is_good": is_good,
            "single_temp": single_temp,
        }
    )
    geds = Table(
        col_dict={
            "energy": energy,
            "rawid": rawid,
            "is_good_channel": is_good_channel,
            "multiplicity": multiplicity,
            "psd": psd,
        }
    )
    # coincident has only a geds flag — no spms → has_spms_coinc=False in pdf.py
    coincident = Table(
        col_dict={"geds": Array(np.array([True, False, True, True, False, True]))}
    )
    evt = Table(col_dict={"geds": geds, "coincident": coincident})
    lh5.write(evt, "evt", path, wo_mode="write_safe")
    _write_cvt_root(path, ["V01", "B02"], [1, 2])


def _make_cvt_file_no_geds(path: Path) -> None:
    """Write a minimal cvt LH5 file with no geds subtable (mimics ``skip_hit: true``)."""
    spms = Array(np.array([False, True, False, True, False, False]))
    coincident = Table(col_dict={"spms": spms})
    # no geds field at the evt level → has_geds=False
    evt = Table(col_dict={"coincident": coincident})
    lh5.write(evt, "evt", path, wo_mode="write_safe")
    # detector_uids is still present (skip_hit produces an empty mapping at evt tier)
    _write_cvt_root(path, [], [])


def test_pdf_script_cli_skip_opt(tmp_path, monkeypatch):
    """With no spms coincidence data (skip_opt), LAr histograms must be absent."""
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file_no_spms(cvt_file)
    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file)

    assert pdf_file.exists()

    inner_keys = lh5.ls(pdf_file, "pdf/")

    # geds-based histograms must all be present
    for name in ("pdf/hit", "pdf/mul", "pdf/mul_psd", "pdf/mul2", "pdf/fail"):
        assert name in inner_keys, f"'{name}' missing from pdf/; got {inner_keys}"

    fail_keys = lh5.ls(pdf_file, "pdf/fail/")
    assert "pdf/fail/psd" in fail_keys, (
        f"'pdf/fail/psd' missing from pdf/fail/; got {fail_keys}"
    )

    # LAr histograms must be absent: no spms coincidence data was present
    for name in ("pdf/mul_lar", "pdf/mul_lar_psd"):
        assert name not in inner_keys, (
            f"'{name}' should be absent when has_spms_coinc=False; got {inner_keys}"
        )
    assert "pdf/fail/lar" not in fail_keys, (
        f"'pdf/fail/lar' should be absent when has_spms_coinc=False; got {fail_keys}"
    )


def test_pdf_script_cli_skip_hit(tmp_path, monkeypatch, caplog):
    """With no geds data (skip_hit), geds histograms absent; LAr histograms present with 0 counts."""
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file_no_geds(cvt_file)
    with caplog.at_level("WARNING"):
        pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file)
    # spurious "matches no detectors" warning is suppressed when has_geds=False
    assert "matches no detectors" not in caplog.text

    assert pdf_file.exists()

    inner_keys = lh5.ls(pdf_file, "pdf/")

    # geds-based histograms must be absent (never initialised when has_geds=False)
    for name in ("pdf/hit", "pdf/mul", "pdf/mul_psd", "pdf/mul2"):
        assert name not in inner_keys, (
            f"'{name}' should be absent when has_geds=False; got {inner_keys}"
        )

    # LAr histograms must be present (initialised from has_spms_coinc=True guard)
    for name in ("pdf/mul_lar", "pdf/mul_lar_psd"):
        assert name in inner_keys, (
            f"'{name}' missing from pdf/; expected when has_spms_coinc=True, got {inner_keys}"
        )

    fail_keys = lh5.ls(pdf_file, "pdf/fail/")
    assert "pdf/fail/lar" in fail_keys, (
        f"'pdf/fail/lar' missing from pdf/fail/; expected when has_spms_coinc=True, got {fail_keys}"
    )

    # fail/psd must be absent (needs geds data)
    assert "pdf/fail/psd" not in fail_keys, (
        f"'pdf/fail/psd' should be absent when has_geds=False; got {fail_keys}"
    )

    def _sum(path):
        return lh5.read_as(path, pdf_file, "hist").sum()

    # all LAr histograms were initialised but never filled → 0 counts
    assert _sum("pdf/mul_lar/all") == 0, (
        "pdf/mul_lar should have 0 counts when has_geds=False"
    )
    assert _sum("pdf/mul_lar_psd/all") == 0, (
        "pdf/mul_lar_psd should have 0 counts when has_geds=False"
    )
    assert _sum("pdf/fail/lar/all") == 0, (
        "pdf/fail/lar should have 0 counts when has_geds=False"
    )


@pytest.mark.needs_remage
def test_pdf_script_cli_with_real_cvt(tmp_path, monkeypatch, legend_cvt_path):
    """Run the pdf script on a real cvt file produced by the full test pipeline."""
    pdf_file = _run_pdf(
        tmp_path,
        monkeypatch,
        legend_cvt_path,
        simid=_SIMID_L1000,
        config_template="simflow-config-l1000.yaml",
    )

    assert pdf_file.exists(), "pdf output file was not created"

    root_keys = lh5.ls(pdf_file)
    assert "pdf" in root_keys
    assert "nr_sim_events" in root_keys

    nr_sim_events = lh5.read("nr_sim_events", pdf_file)
    assert np.issubdtype(type(nr_sim_events.value), np.integer)
    assert nr_sim_events.value > 0

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

    def _sum(path):
        return lh5.read_as(path, pdf_file, "hist").sum()

    # cut hierarchy must be non-increasing
    n_mul = _sum("pdf/mul/all")
    assert _sum("pdf/hit/all") >= n_mul
    assert _sum("pdf/mul_lar/all") <= n_mul
    assert _sum("pdf/mul_psd/all") <= n_mul
    assert _sum("pdf/mul_lar_psd/all") <= _sum("pdf/mul_lar/all")
    assert _sum("pdf/mul_lar_psd/all") <= _sum("pdf/mul_psd/all")


def test_pdf_detector_groups_schema(tmp_path, monkeypatch):
    """With detector_groups configured, output has per-group and 'all' histograms for every cut."""
    meta_dir = _make_metadata_with_pdf_settings(
        tmp_path, {"detector_groups": {"icpc": "V.*", "bege": "B.*"}}
    )

    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)

    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)
    assert pdf_file.exists()

    root_keys = lh5.ls(str(pdf_file))
    assert "pdf" in root_keys
    assert "nr_sim_events" in root_keys
    assert "mul2" not in root_keys  # mul2 is under pdf/, not root

    inner_keys = lh5.ls(str(pdf_file), "pdf/")
    assert "pdf/mul2" in inner_keys  # global, not split by group

    for cut in ("hit", "mul", "mul_lar", "mul_psd", "mul_lar_psd"):
        assert f"pdf/{cut}" in inner_keys, f"pdf/{cut} missing"
        cut_keys = lh5.ls(str(pdf_file), f"pdf/{cut}/")
        for group in ("icpc", "bege", "all"):
            assert f"pdf/{cut}/{group}" in cut_keys, (
                f"pdf/{cut}/{group} missing; got {cut_keys}"
            )

    fail_keys = lh5.ls(str(pdf_file), "pdf/fail/")
    for fail_cut in ("psd", "lar"):
        assert f"pdf/fail/{fail_cut}" in fail_keys, f"pdf/fail/{fail_cut} missing"
        fc_keys = lh5.ls(str(pdf_file), f"pdf/fail/{fail_cut}/")
        for group in ("icpc", "bege", "all"):
            assert f"pdf/fail/{fail_cut}/{group}" in fc_keys, (
                f"pdf/fail/{fail_cut}/{group} missing; got {fc_keys}"
            )


def test_pdf_no_detector_groups_schema(tmp_path, monkeypatch):
    """Without detector_groups, output has only 'all' histograms for every cut."""
    meta_dir = _make_metadata_with_pdf_settings(tmp_path, {})

    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)

    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)
    assert pdf_file.exists()

    inner_keys = lh5.ls(str(pdf_file), "pdf/")
    assert "pdf/mul2" in inner_keys

    for cut in ("hit", "mul", "mul_lar", "mul_psd", "mul_lar_psd"):
        assert f"pdf/{cut}" in inner_keys, f"pdf/{cut} missing"
        cut_keys = lh5.ls(str(pdf_file), f"pdf/{cut}/")
        assert f"pdf/{cut}/all" in cut_keys, f"pdf/{cut}/all missing"
        for group in ("icpc", "bege"):
            assert f"pdf/{cut}/{group}" not in cut_keys, (
                f"pdf/{cut}/{group} should be absent without detector_groups; got {cut_keys}"
            )

    fail_keys = lh5.ls(str(pdf_file), "pdf/fail/")
    for fail_cut in ("psd", "lar"):
        assert f"pdf/fail/{fail_cut}" in fail_keys
        fc_keys = lh5.ls(str(pdf_file), f"pdf/fail/{fail_cut}/")
        assert f"pdf/fail/{fail_cut}/all" in fc_keys
        for group in ("icpc", "bege"):
            assert f"pdf/fail/{fail_cut}/{group}" not in fc_keys, (
                f"pdf/fail/{fail_cut}/{group} should be absent without detector_groups"
            )


def test_pdf_detector_groups_sum_equals_all(tmp_path, monkeypatch):
    """For non-overlapping groups that partition all detectors, sum of groups equals 'all'."""
    meta_dir = _make_metadata_with_pdf_settings(
        tmp_path, {"detector_groups": {"icpc": "V.*", "bege": "B.*"}}
    )

    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)

    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)

    def _counts(path):
        return lh5.read_as(path, str(pdf_file), "hist").values()

    cuts_to_check = [
        "hit",
        "mul",
        "mul_lar",
        "mul_psd",
        "mul_lar_psd",
        "fail/psd",
        "fail/lar",
    ]
    for cut in cuts_to_check:
        icpc = _counts(f"pdf/{cut}/icpc")
        bege = _counts(f"pdf/{cut}/bege")
        all_counts = _counts(f"pdf/{cut}/all")
        np.testing.assert_array_equal(
            icpc + bege,
            all_counts,
            err_msg=f"sum(icpc + bege) != all for cut '{cut}'",
        )
