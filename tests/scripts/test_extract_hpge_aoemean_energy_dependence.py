from __future__ import annotations

import shutil
import sys
from pathlib import Path

import lgdo
import lh5
import numpy as np
import pytest
import yaml
from lgdo import Array, Table, VectorOfVectors

from legendsimflow.metadata import ELECTRON_GUN_ENERGIES_IN_KEV
from legendsimflow.scripts import extract_hpge_aoemean_energy_dependence as script

dummyprod = Path(__file__).parent.parent / "dummyprod"

_RUNID = "l200-p03-r000-phy"
_DET = "V05261B"
# origin of V05261B in the test geometry (m), as written by remage in the
# `detector_origins` struct
_DET_ORIGIN = (0.109, -0.18879354, 0.7951265)
_N_EVENTS = 300


def _write_fake_electron_stp(path: Path, energy_kev: float, seed: int) -> None:
    """Write a remage-like stp file with point-like deposits of `energy_kev` in V05261B.

    Mimics the layout of the electron-gun output: a `stp/{det}` table of events
    with one step each, and the `detector_origins` struct.
    """
    rng = np.random.default_rng(seed)

    # uniform in a cylinder well inside the crystal, in local coordinates (m)
    r = rng.uniform(0.005, 0.025, _N_EVENTS)
    phi = rng.uniform(0, 2 * np.pi, _N_EVENTS)
    z = rng.uniform(0.010, 0.060, _N_EVENTS)

    def _vov(values: np.ndarray) -> VectorOfVectors:
        return VectorOfVectors(
            flattened_data=Array(np.asarray(values, dtype=np.float32)),
            cumulative_length=Array(np.arange(1, _N_EVENTS + 1)),
        )

    table = Table(
        col_dict={
            "evtid": Array(np.arange(_N_EVENTS)),
            "t0": Array(np.zeros(_N_EVENTS), attrs={"units": "ns"}),
            "time": _vov(np.zeros(_N_EVENTS)),
            "edep": _vov(np.full(_N_EVENTS, energy_kev)),
            "xloc": _vov(_DET_ORIGIN[0] + r * np.cos(phi)),
            "yloc": _vov(_DET_ORIGIN[1] + r * np.sin(phi)),
            "zloc": _vov(_DET_ORIGIN[2] + z),
            # far from the surface: the precomputed distance is kept as is
            "dist_to_surf": _vov(np.full(_N_EVENTS, 0.01)),
            "particle": _vov(np.full(_N_EVENTS, 11)),
            "trackid": _vov(np.ones(_N_EVENTS)),
            "parent_trackid": _vov(np.zeros(_N_EVENTS)),
        }
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    lh5.write(table, f"stp/{_DET}", path, wo_mode="write_safe")

    origins = lgdo.Struct(
        {
            _DET: lgdo.Struct(
                {
                    k: lgdo.Scalar(np.float32(v))
                    for k, v in zip(("xloc", "yloc", "zloc"), _DET_ORIGIN, strict=True)
                }
            )
        }
    )
    lh5.write(origins, "detector_origins", path, wo_mode="append")


@pytest.fixture(scope="module")
def electron_stp_files(tmp_path_factory) -> list[Path]:
    """One fake electron-gun stp file per grid energy, named as the par step does."""
    stp_dir = tmp_path_factory.mktemp("electron_stp")
    files = []
    for i, energy in enumerate(ELECTRON_GUN_ENERGIES_IN_KEV):
        f = stp_dir / f"l1000dsg01-electron-gun-{energy}keV-tier_stp.lh5"
        _write_fake_electron_stp(f, energy, seed=i)
        files.append(f)
    return files


def test_group_stp_files_by_energy(electron_stp_files):
    grouped = script.group_stp_files_by_energy(reversed(electron_stp_files))
    assert list(grouped) == sorted(ELECTRON_GUN_ENERGIES_IN_KEV)
    assert all(len(v) == 1 for v in grouped.values())

    with pytest.raises(ValueError):
        script.group_stp_files_by_energy(["/some/birds_nest_K40/file.lh5"])


def test_half_sample_mode():
    rng = np.random.default_rng(1)

    # a Gaussian: mode, median and mean agree
    assert abs(script.half_sample_mode(rng.normal(1.0, 0.01, 20000)) - 1.0) < 2e-3

    # a peak with a heavy low-side tail, as for the electron-gun A/E: the mode
    # tracks the peak while the median is dragged into the tail
    peak = rng.normal(1.0, 0.01, 7000)
    tail = 1.0 - rng.exponential(0.1, 3000)
    sample = np.concatenate([peak, tail])
    mode = script.half_sample_mode(sample)
    assert abs(mode - 1.0) < 5e-3
    assert abs(mode - 1.0) < abs(np.median(sample) - 1.0)

    # degenerate sample sizes
    assert script.half_sample_mode(np.array([3.0])) == 3.0
    assert script.half_sample_mode(np.array([1.0, 3.0])) == 2.0
    assert script.half_sample_mode(np.array([1.0, 1.2, 3.0])) == pytest.approx(1.1)
    with pytest.raises(ValueError):
        script.half_sample_mode(np.array([]))


def test_summarize_aoe():
    rng = np.random.default_rng(1)
    aoe = rng.normal(1.0, 0.01, 10000)
    aoe[:10] = np.nan

    st = script.summarize_aoe(aoe)
    assert st["n_events"] == 9990
    assert abs(st["mode"] - 1.0) < 3e-3
    assert abs(st["median"] - 1.0) < 1e-3
    assert abs(st["std"] - 0.01) < 1e-3
    # the bootstrap error is of the order of the standard error of the mean
    assert 0 < st["mode_err"] < 1e-2
    # MAD-based sigma agrees with the Gaussian sigma
    assert abs(1.4826 * st["mad"] - 0.01) < 1e-3

    empty = script.summarize_aoe(np.array([]))
    assert empty["n_events"] == 0
    assert np.isnan(empty["mode"])


def test_fit_linear_energy_dependence():
    x = np.array([500.0, 1000.0, 1500.0, 2000.0])
    y = 1.0 - 1e-5 * x
    pars, errs = script.fit_linear_energy_dependence(x, y, np.full(4, 1e-4))
    assert abs(pars["a"] + 1e-5) < 1e-8
    assert abs(pars["b"] - 1.0) < 1e-4
    assert errs["a"] > 0
    assert errs["b"] > 0

    # non-finite points are dropped
    pars2, _ = script.fit_linear_energy_dependence(
        [*x, 3000.0], [*y, np.nan], [*np.full(4, 1e-4), 1e-4]
    )
    assert abs(pars2["a"] - pars["a"]) < 1e-9

    with pytest.raises(ValueError):
        script.fit_linear_energy_dependence([500.0], [1.0], [1e-4])


def test_extract_hpge_aoemean_energy_dependence(
    tmp_path,
    monkeypatch,
    electron_stp_files,
    legend_gdml_path,
    legend_dtmap_path,
    legend_currmod_paths,
    l1000_config_factory,
):
    """End-to-end run of the script on point-like deposits.

    Point-like deposits give A/E = 1 (up to the noise smearing) at all energies,
    so the extracted model must be flat and centred at 1. Only V05261B has a
    drift-time map (and hits) in the test inputs: the other modelable detector
    of the run is skipped with a warning.
    """
    # only the single-template PSD: no realistic PSL in the test pars
    config_path = l1000_config_factory(
        tmp_path, settings_by_tier={"hit": {"simulate_psd_with_psl": False}}
    )
    raw = yaml.safe_load(config_path.read_text())
    pars_dir = tmp_path / "pars"
    raw["paths"]["pars"] = str(pars_dir)
    config_path.write_text(yaml.safe_dump(raw))

    dtmap_dir = pars_dir / "hpge/dtmaps"
    dtmap_dir.mkdir(parents=True)
    shutil.copy(legend_dtmap_path, dtmap_dir / f"{_RUNID}-hpge-drift-time-maps.lh5")
    currmod_dir = pars_dir / "hpge/currmod"
    currmod_dir.mkdir(parents=True)
    shutil.copy(legend_currmod_paths[_RUNID], currmod_dir / f"{_RUNID}-model.yaml")

    is_modelable_file = tmp_path / "is_modelable.yaml"
    is_modelable_file.write_text(
        yaml.safe_dump({_RUNID: {"V02160A": True, _DET: True, "B00000A": False}})
    )

    pars_file = tmp_path / "outputs" / f"{_RUNID}-model.yaml"
    stats_file = tmp_path / "outputs" / f"{_RUNID}-electron-gun-stats.yaml"
    plot_file = tmp_path / "plots" / f"{_RUNID}-fit-results.pdf"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "extract-hpge-aoemeanmod",
            "--runid",
            _RUNID,
            "--electron-stp-files",
            *[str(f) for f in electron_stp_files],
            "--geom-file",
            str(legend_gdml_path),
            "--is-modelable-file",
            str(is_modelable_file),
            "--pars-file",
            str(pars_file),
            "--stats-file",
            str(stats_file),
            "--plot-file",
            str(plot_file),
            "--simflow-config",
            str(config_path),
        ],
    )
    script.main()

    assert pars_file.exists()
    assert stats_file.exists()
    assert plot_file.exists()

    # the model file is in the schema consumed by the hit tier
    models = yaml.safe_load(pars_file.read_text())
    assert set(models) == {_DET}
    assert set(models[_DET]) == {"single_template"}
    model = models[_DET]["single_template"]
    assert model["expression"] == "x*a+b"
    assert set(model["pars"]) == {"a", "b"} == set(model["errs"])
    assert np.isfinite(model["pars"]["a"])
    assert np.isfinite(model["pars"]["b"])
    # flat and centred at 1: the slope over the 2.5 MeV lever arm is negligible
    assert abs(model["pars"]["b"] - 1.0) < 0.01
    assert abs(model["pars"]["a"]) * 2500 < 0.01

    stats = yaml.safe_load(stats_file.read_text())
    assert set(stats) == {_DET}
    pts = stats[_DET]["electron_gun"]["single_template"]
    assert pts["energy_in_keV"] == list(ELECTRON_GUN_ENERGIES_IN_KEV)
    assert pts["n_events"] == [_N_EVENTS] * len(ELECTRON_GUN_ENERGIES_IN_KEV)
    assert all(abs(m - 1.0) < 0.01 for m in pts["mode"])
    assert all(e > 0 for e in pts["mode_err"])
    # the noise smearing makes the resolution shrink with energy
    assert pts["std"][0] > pts["std"][-1]

    reso = stats[_DET]["mc_resolution"]["single_template"]
    # reported at the grid energy closest to 2 MeV
    assert reso["energy_in_keV"] == min(
        ELECTRON_GUN_ENERGIES_IN_KEV, key=lambda e: abs(e - 2000)
    )
    assert reso["val"] > 0
