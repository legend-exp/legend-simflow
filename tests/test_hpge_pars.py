from __future__ import annotations

from pathlib import Path

import numpy as np
from lgdo import lh5
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy.stats import norm

from legendsimflow import hpge_pars


def test_fit():
    t = np.linspace(-1000, 2000, 3001)
    A = norm.pdf(t, loc=0, scale=50)

    res = hpge_pars.fit_current_pulse(t, A)

    assert isinstance(res[0], np.ndarray)
    assert isinstance(res[1], np.ndarray)
    assert isinstance(res[2], np.ndarray)


def test_get_index(legend_testdata):
    ref_path = legend_testdata.get_path("lh5/prod-ref-l200/")
    path = ref_path / Path("generated/tier/hit/cal/p03/r001")
    files = [str(p) for p in path.glob("*")]

    # we have to be very generous with the (low stats) test file
    idx, file_idx = hpge_pars.find_best_event_idx(
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


def test_plot():
    fig, ax = hpge_pars.plot_fit_result(
        [1, 2, 3], [0, 10, 20], np.linspace(0, 3, 1000), np.linspace(0, 3, 1000)
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)


def test_get_parameter_dict():
    popt = [100, 10, 60, 0.6, 100, 0.2, 60]

    popt_dict = hpge_pars._curve_fit_popt_to_dict(popt)
    assert len(popt_dict) == 8
