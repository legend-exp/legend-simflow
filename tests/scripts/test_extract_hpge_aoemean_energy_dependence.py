from __future__ import annotations

import shutil
import sys

import numpy as np
import pytest
import yaml

from legendsimflow.metadata import ELECTRON_GUN_ENERGIES_IN_KEV
from legendsimflow.scripts import extract_hpge_aoemean_energy_dependence as script

_RUNID = "l200-p03-r000-phy"
# the two detectors of the mock array that are modelled (they are `on` and have
# a drift-time map)
_DETS = ("V00001A", "V00001B")


def _electron_aoe(rng, n, mu=1.0, sigma=0.01, tail_frac=0.1, tau=0.05):
    """A sample shaped like the electron-gun A/E: a peak with a low-side tail."""
    n_tail = int(n * tail_frac)
    return np.concatenate(
        [
            rng.normal(mu, sigma, n - n_tail),
            rng.normal(mu, sigma, n_tail) - rng.exponential(tau, n_tail),
        ]
    )


def test_fit_aoe_peak():
    rng = np.random.default_rng(1)
    aoe = _electron_aoe(rng, 40000)
    aoe[:10] = np.nan

    stat, curve = script.fit_aoe_peak(aoe)
    assert stat["n_events"] == 39990
    # the Gaussian mean is recovered, and is not dragged into the tail the way
    # the sample mean is
    assert stat["mu"] == pytest.approx(1.0, abs=5e-4)
    assert abs(stat["mu"] - 1.0) < abs(np.nanmean(aoe) - 1.0)
    assert stat["sigma"] == pytest.approx(0.01, rel=0.1)
    assert 0 < stat["mu_err"] < 1e-3
    assert stat["chi2_ndf"] < 3

    # the fitted curve comes back for the plot, on the histogram of the fit
    edges, expected = curve
    assert len(expected) == len(edges) - 1
    assert np.all(np.isfinite(expected))

    # too few events to fit: NaNs, and no curve, so the point is dropped
    stat, curve = script.fit_aoe_peak(np.array([]))
    assert stat["n_events"] == 0
    assert np.isnan(stat["mu"])
    assert curve is None


def test_fit_aoe_vs_energy():
    x = np.array([900.0, 1250.0, 1600.0, 1950.0])
    y = 1.0 - 1e-5 * x

    model = script.fit_aoe_vs_energy(x, y, np.full(4, 1e-4))
    assert model["expression"] == "x*a+b"
    assert abs(model["pars"]["a"] + 1e-5) < 1e-8
    assert abs(model["pars"]["b"] - 1.0) < 1e-4
    assert model["errs"]["a"] > 0
    assert model["errs"]["b"] > 0

    # non-finite points are dropped
    dropped = script.fit_aoe_vs_energy(
        [*x, 2350.0], [*y, np.nan], [*np.full(4, 1e-4), 1e-4]
    )
    assert abs(dropped["pars"]["a"] - model["pars"]["a"]) < 1e-9

    # too few valid points to fit
    assert script.fit_aoe_vs_energy([900.0], [1.0], [1e-4]) is None
    assert script.fit_aoe_vs_energy(x, [np.nan] * 4, np.full(4, 1e-4)) is None
    # a point with no uncertainty carries no information and is dropped
    assert script.fit_aoe_vs_energy(x, y, np.zeros(4)) is None


@pytest.mark.needs_remage
def test_extract_hpge_aoemean_energy_dependence(
    tmp_path,
    monkeypatch,
    legend_electron_gun_paths,
    legend_dtmap_path,
    legend_currmod_paths,
    l200_config_factory,
):
    # only the single-template PSD: no realistic pulse-shape library is built
    # in the tier-script tests
    config_path = l200_config_factory(
        tmp_path, settings_by_tier={"hit": {"simulate_psd_with_psl": False}}
    )
    pars_dir = tmp_path / "generated-l200/pars"

    dtmap_dir = pars_dir / "hpge/dtmaps"
    dtmap_dir.mkdir(parents=True)
    shutil.copy(legend_dtmap_path, dtmap_dir / f"{_RUNID}-hpge-drift-time-maps.lh5")
    currmod_dir = pars_dir / "hpge/currmod"
    currmod_dir.mkdir(parents=True)
    shutil.copy(legend_currmod_paths[_RUNID], currmod_dir / f"{_RUNID}-model.yaml")

    # the detinfo file the cache_modelable_hpges checkpoint writes
    is_modelable_file = tmp_path / "is_modelable.yaml"
    is_modelable_file.write_text(
        yaml.safe_dump({_RUNID: dict.fromkeys(_DETS, True) | {"V00002A": False}})
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
            *[str(f) for f in legend_electron_gun_paths["stp_files"]],
            "--geom-file",
            str(legend_electron_gun_paths["geom_file"]),
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

    assert plot_file.exists()

    # the model file is in the schema consumed by the hit tier
    models = yaml.safe_load(pars_file.read_text())
    assert set(models) == set(_DETS)
    for det in _DETS:
        assert set(models[det]) == {"single_template"}
        model = models[det]["single_template"]
        assert model["expression"] == "x*a+b"
        assert set(model["pars"]) == {"a", "b"} == set(model["errs"])
        # the A/E band of single-site events is centred around 1 and its
        # energy dependence is a small correction over the 1.5 MeV lever arm
        assert model["pars"]["b"] == pytest.approx(1.0, abs=0.1)
        assert abs(model["pars"]["a"]) * 1500 < 0.1

    stats = yaml.safe_load(stats_file.read_text())
    assert set(stats) == set(_DETS)
    for det in _DETS:
        points = stats[det]["single_template"]
        assert points["energy_in_keV"] == list(ELECTRON_GUN_ENERGIES_IN_KEV)
        assert all(n > 0 for n in points["n_events"])
        assert all(abs(m - 1.0) < 0.1 for m in points["mu"])
        assert all(e > 0 for e in points["mu_err"])
        assert all(s > 0 for s in points["sigma"])
        # the full-energy cut leaves a shape the model describes
        assert all(c < 5 for c in points["chi2_ndf"])
        # only the events that stopped the whole electron are used, so fewer
        # than were simulated in the detector
        assert all(n < 40000 for n in points["n_events"])
