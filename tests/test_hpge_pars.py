from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest
from dbetto import AttrsDict
from iminuit import Minuit
from lgdo import lh5
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy.stats import norm

from legendsimflow import hpge_pars


def test_fit():
    t = np.linspace(-1000, 2000, 3001)
    A = norm.pdf(t, loc=0, scale=50)

    res = hpge_pars.fit_currmod(t, A)

    assert isinstance(res[0], np.ndarray)
    assert isinstance(res[1], np.ndarray)
    assert isinstance(res[2], np.ndarray)


def test_fit_gauss():
    # setup the fit
    fitter = hpge_pars.fit_noise_gauss(
        [1, 2, 3], bins=10, fit_range=(0, 5), sigma_range=(2, 20)
    )

    # check fit ran
    assert isinstance(fitter, Minuit)

    # check limits
    assert fitter.values["mu"] >= 0
    assert fitter.values["mu"] <= 5
    assert fitter.values["sigma"] >= 2
    assert fitter.values["sigma"] <= 20

    fig, ax = hpge_pars.plot_gauss_fit([1, 2, 3], fitter, bins=100)

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)


def test_get_index(legend_testdata):
    ref_path = legend_testdata.get_path("lh5/prod-ref-l200/")
    path = ref_path / Path("generated/tier/hit/cal/p03/r001")
    files = [str(p) for p in path.glob("*")]

    # we have to be very generous with the (low stats) test file
    idx, file_idx = hpge_pars.lookup_currmod_fit_data(
        files, "ch1084803/hit", ewin_center=100, ewin_width=20
    )

    assert idx > -1
    assert file_idx > -1

    # now read back in and check

    energy = lh5.read(
        "ch1084803/hit/cuspEmax_ctc_cal", files[file_idx], idx=[idx]
    ).view_as("np")
    AoE = lh5.read("ch1084803/hit/AoE_Classifier", files[file_idx], idx=[idx]).view_as(
        "np"
    )

    assert (energy > 90) & (energy < 110)
    assert abs(AoE) < 1.5


def test_get_waveform(legend_testdata):
    times, wf = hpge_pars.get_current_pulse(
        legend_testdata.get_path(
            "lh5/prod-ref-l200/generated/tier/raw/cal/p03/r001/l200-p03-r001-cal-20230318T012144Z-tier_raw.lh5"
        ),
        "ch1084803/raw",
        idx=0,
        dsp_config=None,
        dsp_output="waveform",
        align=None,
    )

    assert isinstance(times, np.ndarray)
    assert isinstance(wf, np.ndarray)
    assert len(times) == len(wf)


def test_get_waveform_maxima():
    t = np.linspace(-1000, 2000, 3001)

    template = norm.pdf(t, loc=200, scale=50)

    rng = np.random.default_rng()
    noise = rng.uniform(size=(100, 3001), low=0, high=1)

    maxi = hpge_pars.get_waveform_maxima(template, noise)

    assert np.all(maxi > 1)
    assert np.all(maxi < 2)


def test_get_noise_waveforms(legend_testdata):
    raw_file = legend_testdata.get_path(
        "lh5/prod-ref-l200/generated/tier/raw/phy/p03/r001/l200-p03-r001-phy-20230322T160139Z-tier_raw.lh5"
    )
    hit_file = legend_testdata.get_path(
        "lh5/prod-ref-l200/generated/tier/hit/phy/p03/r001/l200-p03-r001-phy-20230322T160139Z-tier_hit.lh5"
    )

    wfs = hpge_pars.get_noise_waveforms(
        [raw_file],
        [hit_file],
        lh5_group="ch1084803/raw",
        energy_var="cuspEmax_ctc_cal",
        dsp_config=None,
        dsp_output="waveform",
    )

    assert wfs is not None
    assert np.shape(wfs)[1] == 1000


def test_get_aoe():
    popt = [100, 10, 60, 0.6, 100, 0.2, 60]

    aoe = hpge_pars.estimate_mean_aoe(popt)
    assert isinstance(aoe, float)


def test_plot():
    fig, ax = hpge_pars.plot_currmod_fit_result(
        [1, 2, 3], [0, 10, 20], np.linspace(0, 3, 1000), np.linspace(0, 3, 1000)
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)


def test_lookup_eres(config, test_l200data):
    meta = hpge_pars.lookup_energy_res_metadata(
        test_l200data / "v2.1.5",
        config.metadata,
        "l200-p03-r000-phy",
        hit_tier_name="pht",
    )

    assert isinstance(meta, AttrsDict)
    for k, v in meta.items():
        assert isinstance(k, str)
        assert "expression" in v
        assert "parameters" in v

    meta = hpge_pars.lookup_energy_res_metadata(
        test_l200data / "v3.0.0",
        config.metadata,
        "l200-p16-r008-ssc",
        hit_tier_name="hit",
    )

    assert isinstance(meta, AttrsDict)
    for k, v in meta.items():
        assert isinstance(k, str)
        assert "expression" in v
        assert "parameters" in v


def test_eres_func():
    f = hpge_pars.build_energy_res_func("FWHMLinear")
    assert isinstance(f, Callable)


def test_build_eres_funcs(config, test_l200data):
    meta = hpge_pars.build_energy_res_func_dict(
        test_l200data / "v2.1.5",
        config.metadata,
        "l200-p03-r000-phy",
        hit_tier_name="pht",
    )

    assert isinstance(meta, dict)
    assert list(meta.keys()) == ["V99000A"]


def test_lookup_aoeres(config, test_l200data):
    meta = hpge_pars.lookup_aoe_res_metadata(
        test_l200data / "v2.1.5",
        config.metadata,
        "l200-p03-r000-phy",
        hit_tier_name="pht",
    )

    assert isinstance(meta, AttrsDict)
    for k, v in meta.items():
        assert isinstance(k, str)
        assert "expression" in v
        assert "pars" in v

    meta = hpge_pars.lookup_aoe_res_metadata(
        test_l200data / "v3.0.0",
        config.metadata,
        "l200-p16-r008-ssc",
        hit_tier_name="hit",
    )

    assert isinstance(meta, AttrsDict)
    for k, v in meta.items():
        assert isinstance(k, str)
        assert "expression" in v
        assert "pars" in v


def test_aoeres_func():
    f = hpge_pars.build_aoe_res_func("SigmaFit")
    assert isinstance(f, Callable)


def test_build_aoeres_funcs(config, test_l200data):
    meta = hpge_pars.build_aoe_res_func_dict(
        test_l200data / "v2.1.5",
        config.metadata,
        "l200-p03-r000-phy",
        hit_tier_name="pht",
    )

    assert isinstance(meta, dict)
    assert list(meta.keys()) == ["V99000A"]
    assert meta["V99000A"](2000) == pytest.approx(0.007, abs=0.001)
    assert meta["V99000A"](2000) == pytest.approx(0.007, abs=0.001)


def test_lookup_psd_cut_vals(config, test_l200data):
    meta = hpge_pars.lookup_psd_cut_values(
        test_l200data / "v2.1.5",
        config.metadata,
        "l200-p03-r000-phy",
        hit_tier_name="pht",
    )

    assert isinstance(meta, AttrsDict)
    for k, v in meta.items():
        assert isinstance(k, str)
        assert "aoe" in v
        assert "low_side" in v.aoe
        assert "high_side" in v.aoe

    meta = hpge_pars.lookup_psd_cut_values(
        test_l200data / "v3.0.0",
        config.metadata,
        "l200-p16-r008-ssc",
        hit_tier_name="hit",
    )

    assert isinstance(meta, AttrsDict)
    for k, v in meta.items():
        assert isinstance(k, str)
        assert "aoe" in v
        assert "low_side" in v.aoe
        assert "high_side" in v.aoe
