from __future__ import annotations

import shutil
import sys
from pathlib import Path

import h5py
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


def _exists(pdf_file: Path, path: str) -> bool:
    with h5py.File(pdf_file, "r") as f:
        return path in f


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
    # - event 1: m1, 1000 keV,     spms=False, PSD valid+simulated+multi-site  → mul1/not_aoe_st
    # - event 2: m1, 1500 keV,     spms=False, PSD not valid                   → cut (no observable)
    # - event 3: m1, 2000 keV,     spms=True,  PSD valid+simulated+single-site → mul1/not_lar
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
            "quality": Table(col_dict={"is_good_channel": is_good_channel}),
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

    for name in (
        "pdf/hit",
        "pdf/mul1",
        "pdf/mul1/lar",
        "pdf/mul1/not_lar",
        "pdf/mul1/aoe_st",
        "pdf/mul1/not_aoe_st",
        "pdf/mul1/lar/aoe_st",
        "pdf/mul2",
    ):
        assert _exists(pdf_file, name), name

    def _sum(path):
        return lh5.read_as(path, pdf_file, "hist").sum()

    # 7 per-detector deposits: 5 m1 + 2 m2
    assert _sum("pdf/hit/all") == 7

    # 5 m1 events
    assert _sum("pdf/mul1/all") == 5

    # event 3 (spms=True) is vetoed: 4 survive
    assert _sum("pdf/mul1/lar/all") == 4

    # events 0 and 3 pass PSD (valid+simulated+SS); events 1 (MS), 2 (invalid), 5 (not simulated) cut
    assert _sum("pdf/mul1/aoe_st/all") == 2

    # event 3 also fails LAr, so only event 0 survives both
    assert _sum("pdf/mul1/lar/aoe_st/all") == 1

    # 1 m2 event
    assert _sum("pdf/mul2") == 1

    # event 3 fails LAr
    assert _sum("pdf/mul1/not_lar/all") == 1

    # event 1 (valid+simulated PSD, multi-site) fails PSD; event 5 (not simulated) excluded
    assert _sum("pdf/mul1/not_aoe_st/all") == 1


def test_pdf_pulse_lib_cuts(tmp_path, monkeypatch):
    # all events m1; columns: has_aoe, is_single_site, is_high_aoe, is_bb_like
    # - event 0: 500 keV,  T, T, F, T, spms=False -> passes all cuts
    # - event 1: 1000 keV, T, F, F, F, spms=False -> fails low side
    # - event 2: 1500 keV, T, T, T, F, spms=False -> fails high side
    # - event 3: 2000 keV, F, F, F, F, spms=False -> no A/E, excluded
    # - event 4: 2500 keV, T, T, F, T, spms=True  -> passes all, LAr vetoed
    def _vov(values):
        return VectorOfVectors(data=[[v] for v in values])

    pulse_lib = Table(
        col_dict={
            "has_aoe": _vov([True, True, True, False, True]),
            "is_single_site": _vov([True, False, True, False, True]),
            "is_high_aoe": _vov([False, False, True, False, False]),
            "is_bb_like": _vov([True, False, False, False, True]),
        }
    )
    psd = Table(col_dict={"is_good": _vov([True] * 5), "pulse_lib": pulse_lib})
    geds = Table(
        col_dict={
            "energy": _vov([500.0, 1000.0, 1500.0, 2000.0, 2500.0]),
            "rawid": _vov([1] * 5),
            "quality": Table(col_dict={"is_good_channel": _vov([True] * 5)}),
            "multiplicity": Array(np.ones(5, dtype=np.int32)),
            "psd": psd,
        }
    )
    spms = Array(np.array([False, False, False, False, True]))
    evt = Table(col_dict={"geds": geds, "coincident": Table(col_dict={"spms": spms})})

    cvt_file = tmp_path / "cvt.lh5"
    lh5.write(evt, "evt", cvt_file, wo_mode="write_safe")
    _write_cvt_root(cvt_file, ["V01"], [1])
    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file)

    # no single_temp model in the input: no single-template histograms
    assert not _exists(pdf_file, "pdf/mul1/aoe_st")
    assert not _exists(pdf_file, "pdf/mul1/not_aoe_st")

    def _sum(path):
        return lh5.read_as(path, pdf_file, "hist").sum()

    expected = {
        "aoe_psl_low": (3, 2, 1),
        "aoe_psl_high": (3, 2, 1),
        "aoe_psl": (2, 1, 2),
    }
    for cut, (n_pass, n_pass_lar, n_fail) in expected.items():
        assert _sum(f"pdf/mul1/{cut}/all") == n_pass
        assert _sum(f"pdf/mul1/lar/{cut}/all") == n_pass_lar
        assert _sum(f"pdf/mul1/not_{cut}/all") == n_fail


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
            "quality": Table(col_dict={"is_good_channel": is_good_channel}),
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

    # geds-based histograms must all be present
    for name in (
        "pdf/hit",
        "pdf/mul1",
        "pdf/mul1/aoe_st",
        "pdf/mul1/not_aoe_st",
        "pdf/mul2",
    ):
        assert _exists(pdf_file, name), name

    # LAr histograms must be absent: no spms coincidence data was present
    for name in ("pdf/mul1/lar", "pdf/mul1/not_lar"):
        assert not _exists(pdf_file, name), name


def test_pdf_script_cli_skip_hit(tmp_path, monkeypatch, caplog):
    """With no geds data (skip_hit), no histogram can be filled: none must be written."""
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file_no_geds(cvt_file)
    with caplog.at_level("WARNING"):
        pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file)
    # spurious "matches no detectors" warning is suppressed when has_geds=False
    assert "matches no detectors" not in caplog.text

    assert pdf_file.exists()

    # no histogram is filled without HPGe data, including the LAr-veto ones
    for name in ("pdf/hit", "pdf/mul1", "pdf/mul2"):
        assert not _exists(pdf_file, name), name


@pytest.mark.needs_remage
def test_pdf_script_cli_with_real_cvt(tmp_path, monkeypatch, legend_cvt_path):
    """Run the pdf script on a real cvt file produced by the full test pipeline."""
    pdf_file = _run_pdf(
        tmp_path,
        monkeypatch,
        legend_cvt_path,
        simid=_SIMID_L1000,
        config_template="simflow-config-l200.yaml",
    )

    assert pdf_file.exists(), "pdf output file was not created"

    root_keys = lh5.ls(pdf_file)
    assert "pdf" in root_keys
    assert "nr_sim_events" in root_keys

    nr_sim_events = lh5.read("nr_sim_events", pdf_file)
    assert np.issubdtype(type(nr_sim_events.value), np.integer)
    assert nr_sim_events.value > 0

    for name in (
        "pdf/hit",
        "pdf/mul1",
        "pdf/mul1/lar",
        "pdf/mul1/aoe_st",
        "pdf/mul1/lar/aoe_st",
        "pdf/mul2",
    ):
        assert _exists(pdf_file, name), name

    def _sum(path):
        return lh5.read_as(path, pdf_file, "hist").sum()

    # cut hierarchy must be non-increasing
    n_mul = _sum("pdf/mul1/all")
    assert _sum("pdf/hit/all") >= n_mul
    assert _sum("pdf/mul1/lar/all") <= n_mul
    assert _sum("pdf/mul1/aoe_st/all") <= n_mul
    assert _sum("pdf/mul1/lar/aoe_st/all") <= _sum("pdf/mul1/lar/all")
    assert _sum("pdf/mul1/lar/aoe_st/all") <= _sum("pdf/mul1/aoe_st/all")


_ALL_CUTS = (
    "hit",
    "mul1",
    "mul1/lar",
    "mul1/not_lar",
    "mul1/aoe_st",
    "mul1/not_aoe_st",
    "mul1/lar/aoe_st",
)


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

    assert _exists(pdf_file, "pdf/mul2")  # global, not split by group

    for cut in _ALL_CUTS:
        for group in ("icpc", "bege", "all"):
            assert _exists(pdf_file, f"pdf/{cut}/{group}"), f"pdf/{cut}/{group}"


def test_pdf_no_detector_groups_schema(tmp_path, monkeypatch):
    """Without detector_groups, output has only 'all' histograms for every cut."""
    meta_dir = _make_metadata_with_pdf_settings(tmp_path, {})

    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)

    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)
    assert pdf_file.exists()

    assert _exists(pdf_file, "pdf/mul2")

    for cut in _ALL_CUTS:
        assert _exists(pdf_file, f"pdf/{cut}/all"), f"pdf/{cut}/all"
        for group in ("icpc", "bege"):
            assert not _exists(pdf_file, f"pdf/{cut}/{group}"), f"pdf/{cut}/{group}"


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

    for cut in _ALL_CUTS:
        icpc = _counts(f"pdf/{cut}/icpc")
        bege = _counts(f"pdf/{cut}/bege")
        all_counts = _counts(f"pdf/{cut}/all")
        np.testing.assert_array_equal(
            icpc + bege,
            all_counts,
            err_msg=f"sum(icpc + bege) != all for cut '{cut}'",
        )


@pytest.mark.parametrize(("all_regex", "warns"), [(".*", False), ("V.*", True)])
def test_pdf_explicit_all_group(tmp_path, monkeypatch, caplog, all_regex, warns):
    meta_dir = _make_metadata_with_pdf_settings(
        tmp_path, {"detector_groups": {"all": all_regex, "icpc": "V.*"}}
    )
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)

    with caplog.at_level("WARNING"):
        _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)
    assert ("overridden by the implicit all-detector group" in caplog.text) == warns


def test_pdf_group_name_reserved(tmp_path, monkeypatch):
    meta_dir = _make_metadata_with_pdf_settings(
        tmp_path, {"detector_groups": {"not_lar": "V.*"}}
    )
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)

    with pytest.raises(ValueError, match="reserved"):
        _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)
