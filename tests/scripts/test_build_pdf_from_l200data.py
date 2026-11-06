from __future__ import annotations

import shutil
import sys
from pathlib import Path

import lh5
import numpy as np
import pytest
import yaml
from lgdo import Array, Table, VectorOfVectors

from legendsimflow.exceptions import SimflowConfigError
from legendsimflow.scripts import build_pdf_from_l200data

dummyprod = Path(__file__).parent.parent / "dummyprod"
test_l200data = Path(__file__).parent.parent / "l200data" / "v3.0.0"

# rawids of V99000A and B99000A in the p02 channel map of dummyprod
V, B = 1234568, 1234569


def _vov(values):
    return VectorOfVectors(data=values)


def _hits(values):
    """One hit per event, two in events 4 and 8, both with the same flag."""
    return _vov([[v, v] if i in (4, 8) else [v] for i, v in enumerate(values)])


def _make_pet_file(path: Path) -> None:
    # - event 0: m1, 500 keV,   spms=False, single-site         -> not_spms_coincident/is_ss_aoe
    # - event 1: m1, 1000 keV,  spms=False, multi-site          -> not_spms_coincident/is_low_aoe
    # - event 2: m1, 1500 keV,  spms=False, A/E not good        -> not_spms_coincident/is_invalid_aoe
    # - event 3: m1, 2000 keV,  spms=True,  above high-side cut -> is_spms_coincident/is_high_aoe (pulse_lib)
    # - event 4: m2, 200+300 keV
    # - event 5: m1, 2500 keV,  forced trigger                  -> dropped
    # - event 6: m1, 600 keV in B99000A, not modeled            -> not_spms_coincident/is_invalid_aoe
    # - event 7: m1, 700 keV,   pulser                          -> dropped
    # - event 8: m1, 800 keV + 10 keV below threshold in B      -> not_spms_coincident/is_ss_aoe
    # - event 9: m1, 900 keV,   not bb-like                     -> dropped
    T, F = True, False
    energy = [500.0, 1000.0, 1500.0, 2000.0, 200.0, 2500.0, 600.0, 700.0, 800.0, 900.0]
    energy = _vov(
        [
            [e, 300.0] if i == 4 else [e, 10.0] if i == 8 else [e]
            for i, e in enumerate(energy)
        ]
    )
    geds = Table(
        col_dict={
            "energy": energy,
            "rawid": _vov([[V], [V], [V], [V], [V, B], [V], [B], [V], [V, B], [V]]),
            "multiplicity": Array(np.array([1, 1, 1, 1, 2, 1, 1, 1, 1, 1])),
            "quality": Table(
                col_dict={
                    "is_good_channel": _hits([T] * 10),
                    "is_bb_like": Array(np.array([T] * 9 + [F])),
                }
            ),
            "psd": Table(
                col_dict={
                    # all PSD classifiers: not read, the A/E validity is low_aoe.is_good
                    "is_good": _hits([F] * 10),
                    "low_aoe": Table(
                        col_dict={
                            "is_good": _hits([T, T, F, T, T, T, T, T, T, T]),
                            "is_single_site": _hits([T, F, F, T, T, T, T, T, T, T]),
                        }
                    ),
                    "high_aoe": Table(
                        col_dict={"is_bulk": _hits([T, T, T, F, T, T, T, T, T, T])}
                    ),
                }
            ),
        }
    )
    forced = np.array([F, F, F, F, F, T, F, F, F, F])
    puls = np.array([F, F, F, F, F, F, F, T, F, F])
    # event 3 fails the LAr veto with 5 pe summed; the thresholds are 4
    # channels and 4 pe
    energy_sum = np.full(10, 0.5)
    energy_sum[3] = 5
    spms = Table(
        col_dict={
            "multiplicity": Array(np.full(10, 2)),
            "energy_sum": Array(energy_sum),
        }
    )
    evt = Table(
        col_dict={
            "trigger": Table(
                col_dict={
                    "period": Array(np.full(10, 2)),
                    "run": Array(np.zeros(10, dtype=int)),
                    "is_forced": Array(forced),
                }
            ),
            "coincident": Table(
                col_dict={
                    "puls": Array(puls),
                    "muon": Array(np.zeros(10, dtype=bool)),
                    "muon_offline": Array(np.zeros(10, dtype=bool)),
                }
            ),
            "spms": spms,
            "geds": geds,
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    lh5.write(evt, "evt", path, wo_mode="write_safe")


def _run(
    tmp_path: Path,
    monkeypatch,
    pulse_lib: bool = False,
    lookup: bool = False,
    defaults: bool = False,
    runid: str = "l200-p02-r000-phy",
    evt_settings: dict | None = None,
):
    """Run the script on a pet file of run l200-p02-r000-phy of the dummy production."""
    meta_dir = tmp_path / "inputs"
    shutil.copytree(dummyprod / "inputs", meta_dir)
    overrides = {"evt": evt_settings or {}}
    if pulse_lib:
        overrides["hit"] = {"simulate_psd": False, "simulate_psd_with_psl": True}
    for tier, values in overrides.items():
        f = meta_dir / f"simprod/config/tier/{tier}/legend/settings.yaml"
        settings = yaml.safe_load(f.read_text())
        settings.update(values)
        f.write_text(yaml.safe_dump(settings))

    l200data = tmp_path / "l200data"
    l200data.mkdir()
    shutil.copy(test_l200data / "dataflow-config.yaml", l200data)
    pet_file = (
        l200data
        / "generated/tier/pet/phy/p02/r000/l200-p02-r000-phy-20220102T000000Z-tier_pet.lh5"
    )
    _make_pet_file(pet_file)

    config_path = tmp_path / "simflow-config.yaml"
    raw_config = yaml.safe_load((dummyprod / "simflow-config.yaml").read_text())
    raw_config["paths"]["metadata"] = str(meta_dir)
    raw_config["paths"]["l200data"] = str(l200data)
    raw_config["paths"]["pars"] = str(tmp_path / "pars")
    raw_config["paths"]["tier"]["pdf"] = str(tmp_path / "pdf")
    config_path.write_text(yaml.safe_dump(raw_config))

    argv = ["tier-pdf-l200data", "--runid", runid, "--simflow-config", str(config_path)]
    pdf_file = tmp_path / "pdf.lh5"
    if defaults:
        # the modeled detectors from the cache of the par step, the output in
        # the pdf tier of the production
        cache = tmp_path / "pars/detinfo/is_modelable.yaml"
        cache.parent.mkdir(parents=True)
        cache.write_text(
            yaml.safe_dump({"l200-p02-r000-phy": {"V99000A": True, "B99000A": False}})
        )
        pdf_file = tmp_path / "pdf/l200data" / f"{runid}-tier_pdf.lh5"
    else:
        argv += ["--modeled-detectors", "V99000A", "--pdf-file", str(pdf_file)]
    if not lookup:
        argv += ["--evt-files", str(pet_file)]
    monkeypatch.setattr(sys, "argv", argv)
    build_pdf_from_l200data.main()
    return pdf_file


def _sum(pdf_file: Path, path: str) -> float:
    return lh5.read_as(path, pdf_file, "hist").sum()


@pytest.mark.parametrize(
    ("lookup", "defaults"), [(False, False), (True, False), (True, True)]
)
def test_build_pdf_from_l200data(tmp_path, monkeypatch, lookup, defaults):
    pdf_file = _run(tmp_path, monkeypatch, lookup=lookup, defaults=defaults)

    assert lh5.read("livetime_in_s", pdf_file).value == 10000
    expected = {
        "not_spms_coincident/is_ss_aoe": 2,
        "not_spms_coincident/is_low_aoe": 1,
        "not_spms_coincident/is_invalid_aoe": 2,
        # single template: no high-side cut, event 3 is single-site
        "is_spms_coincident/is_ss_aoe": 1,
        "is_spms_coincident/is_low_aoe": 0,
        "is_spms_coincident/is_invalid_aoe": 0,
    }
    for path, n in expected.items():
        assert _sum(pdf_file, f"pdf/mul1/{path}/all") == n, path
    assert _sum(pdf_file, "pdf/mul2/all") == 1
    assert (
        lh5.read_as("pdf/mul1/not_spms_coincident/is_ss_aoe/all", pdf_file, "hist")[
            800j
        ]
        == 1
    )


def test_build_pdf_from_l200data_pulse_lib(tmp_path, monkeypatch):
    pdf_file = _run(tmp_path, monkeypatch, pulse_lib=True)

    assert _sum(pdf_file, "pdf/mul1/is_spms_coincident/is_high_aoe/all") == 1
    assert _sum(pdf_file, "pdf/mul1/is_spms_coincident/is_ss_aoe/all") == 0
    assert _sum(pdf_file, "pdf/mul1/not_spms_coincident/is_high_aoe/all") == 0
    assert _sum(pdf_file, "pdf/mul1/not_spms_coincident/is_ss_aoe/all") == 2


def test_build_pdf_from_l200data_run_not_in_cache(tmp_path, monkeypatch):
    with pytest.raises(SimflowConfigError, match="not in the modelable-HPGe cache"):
        _run(tmp_path, monkeypatch, defaults=True, runid="l200-p02-r001-phy")


def test_build_pdf_from_l200data_threshold_below_data(tmp_path, monkeypatch):
    with pytest.raises(SimflowConfigError, match="no hits below 25 keV"):
        _run(tmp_path, monkeypatch, evt_settings={"geds_energy_thr_kev": 20})
