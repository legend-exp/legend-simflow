from __future__ import annotations

import shutil
import sys
from pathlib import Path

import awkward as ak
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


def _sum(pdf_file: Path, path: str) -> float:
    return lh5.read_as(path, pdf_file, "hist").sum()


def _make_metadata(tmp_path: Path, **tier_settings: dict) -> Path:
    """Copy dummyprod inputs to tmp_path and patch tier/<tier>/legend/settings.yaml."""
    meta_dir = tmp_path / "inputs"
    shutil.copytree(dummyprod / "inputs", meta_dir)
    for tier, extra in tier_settings.items():
        path = (
            meta_dir / "simprod" / "config" / "tier" / tier / "legend" / "settings.yaml"
        )
        base = yaml.safe_load(path.read_text()) if path.exists() else {}
        base.update(extra)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(base))
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


def _make_cvt_file(path: Path, with_psd: bool = True) -> None:
    """Write a minimal cvt LH5 file with the fields expected by the pdf script."""
    # 6 events:
    # - event 0: m1, 500 keV,      spms=False, PSD valid+simulated+single-site → not_spms_coincident/is_ss_aoe
    # - event 1: m1, 1000 keV,     spms=False, PSD valid+simulated+multi-site  → not_spms_coincident/is_low_aoe
    # - event 2: m1, 1500 keV,     spms=False, PSD not valid                   → not_spms_coincident/is_invalid_aoe
    # - event 3: m1, 2000 keV,     spms=True,  PSD valid+simulated+single-site → is_spms_coincident/is_ss_aoe
    # - event 4: m2, 200+300 keV,  spms=False
    # - event 5: m1, 2500 keV,     spms=False, PSD valid+not simulated         → not_spms_coincident/is_invalid_aoe
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
        }
    )
    if with_psd:
        geds.add_field("psd", psd)
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

    expected = {
        "not_spms_coincident/is_ss_aoe": 1,
        "not_spms_coincident/is_low_aoe": 1,
        "not_spms_coincident/is_invalid_aoe": 2,
        "is_spms_coincident/is_ss_aoe": 1,
        "is_spms_coincident/is_low_aoe": 0,
        "is_spms_coincident/is_invalid_aoe": 0,
    }
    for path, n in expected.items():
        assert _sum(pdf_file, f"pdf/mul1/{path}/all") == n, path

    # single template: no high-side cut
    assert not _exists(pdf_file, "pdf/mul1/not_spms_coincident/is_high_aoe")
    assert not _exists(pdf_file, "pdf/mul1/all")
    assert not _exists(pdf_file, "pdf/hit")

    # 1 m2 event
    assert _sum(pdf_file, "pdf/mul2/all") == 1


def test_pdf_pulse_lib_cuts(tmp_path, monkeypatch):
    # all events m1; columns: has_aoe, is_single_site, is_high_aoe
    # - event 0: 500 keV,  T, T, F, spms=False -> is_ss_aoe
    # - event 1: 1000 keV, T, F, F, spms=False -> is_low_aoe
    # - event 2: 1500 keV, T, T, T, spms=False -> is_high_aoe
    # - event 3: 2000 keV, F, F, F, spms=False -> is_invalid_aoe
    # - event 4: 2500 keV, T, T, F, spms=True  -> is_ss_aoe
    def _vov(values):
        return VectorOfVectors(data=[[v] for v in values])

    pulse_lib = Table(
        col_dict={
            "has_aoe": _vov([True, True, True, False, True]),
            "is_single_site": _vov([True, False, True, False, True]),
            "is_high_aoe": _vov([False, False, True, False, False]),
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
    meta_dir = _make_metadata(
        tmp_path, hit={"simulate_psd": False, "simulate_psd_with_psl": True}
    )
    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)

    # (not_spms_coincident, is_spms_coincident)
    expected = {
        "is_high_aoe": (1, 0),
        "is_ss_aoe": (1, 1),
        "is_low_aoe": (1, 0),
        "is_invalid_aoe": (1, 0),
    }
    for aoe, (n_no_light, n_light) in expected.items():
        assert (
            _sum(pdf_file, f"pdf/mul1/not_spms_coincident/{aoe}/all") == n_no_light
        ), aoe
        assert _sum(pdf_file, f"pdf/mul1/is_spms_coincident/{aoe}/all") == n_light, aoe


def test_pdf_no_psd(tmp_path, monkeypatch):
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file, with_psd=False)
    meta_dir = _make_metadata(tmp_path, hit={"simulate_psd": False})
    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)

    # without PSD the A/E level is left out: events 0, 1, 2, 5 pass the LAr
    # veto, event 3 does not
    assert _sum(pdf_file, "pdf/mul1/not_spms_coincident/all") == 4
    assert _sum(pdf_file, "pdf/mul1/is_spms_coincident/all") == 1
    assert not _exists(pdf_file, "pdf/mul1/not_spms_coincident/is_ss_aoe")


def _make_cvt_file_surface(path: Path) -> None:
    """Write a cvt file with surface steps and hits below the 25 keV threshold."""
    # detector V99000A has a FCCD of 0.75 mm; a "-" hit is below threshold and
    # in an AC detector, it must be ignored
    # - event 0: m1, 1000 keV, 5 keV in the FCCD             -> not_fccd_tagged
    # - event 1: m1, 1000 keV, 15 keV in the FCCD            -> is_fccd_tagged
    # - event 2: m1, 800 keV + 10 keV(-) with 50 keV in FCCD -> not_fccd_tagged
    # - event 3: m2, 300 + 400 keV + 5 keV(-)
    geds = Table(
        col_dict={
            "energy": VectorOfVectors(
                data=[[1000.0], [1000.0], [800.0, 10.0], [300.0, 400.0, 5.0]]
            ),
            "rawid": VectorOfVectors(data=[[1], [1], [1, 2], [1, 2, 2]]),
            "quality": Table(
                col_dict={
                    "is_good_channel": VectorOfVectors(
                        data=[[True], [True], [True, False], [True, True, False]]
                    )
                }
            ),
            "multiplicity": Array(np.array([1, 1, 1, 2], dtype=np.int32)),
            "surface": Table(
                col_dict={
                    "edep": VectorOfVectors(
                        ak.Array(
                            [[[5.0, 20.0]], [[15.0]], [[1.0], [50.0]], [[], [], []]]
                        )
                    ),
                    "distance_to_nplus": VectorOfVectors(
                        ak.Array([[[0.5, 2.0]], [[0.3]], [[1.0], [0.1]], [[], [], []]])
                    ),
                }
            ),
        }
    )
    lh5.write(Table(col_dict={"geds": geds}), "evt", path, wo_mode="write_safe")
    _write_cvt_root(path, ["V99000A", "B99000A"], [1, 2])


def test_pdf_energy_threshold(tmp_path, monkeypatch):
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file_surface(cvt_file)
    meta_dir = _make_metadata(
        tmp_path,
        evt={"skip_opt": True},
        hit={"simulate_psd": False, "store_surface_steps": True},
    )
    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)

    def _hist(path):
        return lh5.read_as(path, pdf_file, "hist")

    # no SiPM or PSD data: no category levels
    mul1 = _hist("pdf/mul1/all")
    assert mul1.sum() == 3
    assert mul1[10j] == 0
    assert mul1[800j] == 1

    mul2 = _hist("pdf/mul2/all")
    assert mul2.sum() == 1
    assert mul2[300j, 400j] == 1


def test_pdf_fccd_tag(tmp_path, monkeypatch):
    meta_dir = _make_metadata(
        tmp_path,
        evt={"skip_opt": True},
        hit={"simulate_psd": False, "store_surface_steps": True},
        pdf={"fccd_tag_energy_thr_kev": 10},
    )
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file_surface(cvt_file)
    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)

    assert _sum(pdf_file, "pdf/mul1/is_fccd_tagged/all") == 1
    assert _sum(pdf_file, "pdf/mul1/not_fccd_tagged/all") == 2


def test_pdf_fccd_tag_needs_surface_steps(tmp_path, monkeypatch):
    meta_dir = _make_metadata(tmp_path, pdf={"fccd_tag_energy_thr_kev": 10})
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)

    with pytest.raises(ValueError, match="store_surface_steps"):
        _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)


def _make_cvt_file_no_geds(path: Path) -> None:
    """Write a minimal cvt LH5 file with no geds subtable (mimics ``skip_hit: true``)."""
    spms = Array(np.array([False, True, False, True, False, False]))
    coincident = Table(col_dict={"spms": spms})
    evt = Table(col_dict={"coincident": coincident})
    lh5.write(evt, "evt", path, wo_mode="write_safe")
    # detector_uids is still present (skip_hit produces an empty mapping at evt tier)
    _write_cvt_root(path, [], [])


def test_pdf_script_cli_skip_opt(tmp_path, monkeypatch):
    """With no spms coincidence data (skip_opt), LAr histograms must be absent."""
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)
    meta_dir = _make_metadata(tmp_path, evt={"skip_opt": True})
    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)

    assert pdf_file.exists()

    # the LAr level is omitted
    for aoe in ("is_ss_aoe", "is_low_aoe", "is_invalid_aoe"):
        assert _exists(pdf_file, f"pdf/mul1/{aoe}/all"), aoe
    assert _exists(pdf_file, "pdf/mul2/all")
    for name in ("pdf/mul1/is_spms_coincident", "pdf/mul1/not_spms_coincident"):
        assert not _exists(pdf_file, name), name


def test_pdf_script_cli_skip_hit(tmp_path, monkeypatch, caplog):
    """With no geds data (skip_hit), no histogram can be filled: none must be written."""
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file_no_geds(cvt_file)
    with caplog.at_level("WARNING"):
        pdf_file = _run_pdf(
            tmp_path,
            monkeypatch,
            cvt_file,
            meta_dir=_make_metadata(tmp_path, evt={"skip_hit": True}),
        )
    # no detector groups are resolved without HPGe data
    assert "matches no detectors" not in caplog.text

    assert pdf_file.exists()

    # no histogram is filled without HPGe data, including the LAr-veto ones
    for name in ("pdf/mul1", "pdf/mul2"):
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
        "pdf/mul1/is_spms_coincident/is_ss_aoe/all",
        "pdf/mul1/not_spms_coincident/is_ss_aoe/all",
        "pdf/mul2/all",
    ):
        assert _exists(pdf_file, name), name


_ALL_CATEGORIES = tuple(
    f"mul1/{lar}/{aoe}"
    for lar in ("is_spms_coincident", "not_spms_coincident")
    for aoe in ("is_ss_aoe", "is_low_aoe", "is_invalid_aoe")
)


def test_pdf_detector_groups_schema(tmp_path, monkeypatch):
    """With detector_groups configured, output has per-group and 'all' histograms for every category."""
    meta_dir = _make_metadata(
        tmp_path, pdf={"detector_groups": {"icpc": "V.*", "bege": "B.*"}}
    )

    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)

    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)
    assert pdf_file.exists()

    root_keys = lh5.ls(str(pdf_file))
    assert "pdf" in root_keys
    assert "nr_sim_events" in root_keys
    assert "mul2" not in root_keys  # mul2 is under pdf/, not root

    # not split by group
    assert _exists(pdf_file, "pdf/mul2/all")
    assert not _exists(pdf_file, "pdf/mul2/icpc")

    for path in _ALL_CATEGORIES:
        for group in ("icpc", "bege", "all"):
            assert _exists(pdf_file, f"pdf/{path}/{group}"), f"pdf/{path}/{group}"


def test_pdf_no_detector_groups_schema(tmp_path, monkeypatch):
    """Without detector_groups, output has only 'all' histograms for every category."""
    meta_dir = _make_metadata(tmp_path, pdf={})

    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)

    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)
    assert pdf_file.exists()

    assert _exists(pdf_file, "pdf/mul2/all")

    for path in _ALL_CATEGORIES:
        assert _exists(pdf_file, f"pdf/{path}/all"), f"pdf/{path}/all"
        for group in ("icpc", "bege"):
            assert not _exists(pdf_file, f"pdf/{path}/{group}"), f"pdf/{path}/{group}"


def test_pdf_detector_groups_sum_equals_all(tmp_path, monkeypatch):
    """For non-overlapping groups that partition all detectors, sum of groups equals 'all'."""
    meta_dir = _make_metadata(
        tmp_path, pdf={"detector_groups": {"icpc": "V.*", "bege": "B.*"}}
    )

    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)

    pdf_file = _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)

    def _counts(path):
        return lh5.read_as(path, str(pdf_file), "hist").values()

    for path in _ALL_CATEGORIES:
        icpc = _counts(f"pdf/{path}/icpc")
        bege = _counts(f"pdf/{path}/bege")
        all_counts = _counts(f"pdf/{path}/all")
        np.testing.assert_array_equal(
            icpc + bege,
            all_counts,
            err_msg=f"sum(icpc + bege) != all for '{path}'",
        )


@pytest.mark.parametrize(("all_regex", "warns"), [(".*", False), ("V.*", True)])
def test_pdf_explicit_all_group(tmp_path, monkeypatch, caplog, all_regex, warns):
    meta_dir = _make_metadata(
        tmp_path, pdf={"detector_groups": {"all": all_regex, "icpc": "V.*"}}
    )
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_file(cvt_file)

    with caplog.at_level("WARNING"):
        _run_pdf(tmp_path, monkeypatch, cvt_file, meta_dir=meta_dir)
    assert ("overridden by the implicit all-detector group" in caplog.text) == warns
