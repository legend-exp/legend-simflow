from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from lgdo import lh5
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy.stats import norm

from legendsimflow.hpge_pars import (
    fit_waveform,
    get_index,
    get_parameter_dict,
    get_raw_file,
    get_waveform,
    plot_fit,
)


def test_fit():
    t = np.linspace(-1000, 2000, 3001)
    A = norm.pdf(t, loc=0, scale=50)

    res = fit_waveform(t, A)

    assert isinstance(res[0], np.ndarray)
    assert isinstance(res[1], np.ndarray)
    assert isinstance(res[2], np.ndarray)


def test_get_index(legend_testdata):
    ref_path = legend_testdata.get_path("lh5/prod-ref-l200/")
    path = ref_path / Path("generated/tier/hit/cal/p03/r001")
    files = [str(p) for p in path.glob("*")]

    # we have to be very generous with the (low stats) test file
    idx, file_idx = get_index(files, "ch1084803/hit", peak=100, peak_range=10)

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


@pytest.fixture(scope="session")
def test_get_rawfile(legend_testdata):
    ref_path = legend_testdata.get_path("lh5/prod-ref-l200/")
    path = ref_path / Path("generated/tier/hit/")
    raw_path = ref_path / Path("generated/tier/raw/")

    hit_file = next(p for p in path.glob("cal/p03/r001/*"))

    file = get_raw_file(hit_file, path, raw_path, hit_name="hit")

    # check file is lh5
    assert file.suffix == ".lh5"

    return file


def test_get_waveform(test_get_rawfile):
    times, wf = get_waveform(
        test_get_rawfile,
        "ch1084803/raw",
        idx=0,
        dsp_config=None,
        align=None,
        current="waveform",
    )

    assert isinstance(times, np.ndarray)
    assert isinstance(wf, np.ndarray)
    assert len(times) == len(wf)


def test_plot():
    fig, ax = plot_fit(
        [1, 2, 3], [0, 10, 20], np.linspace(0, 3, 1000), np.linspace(0, 3, 1000)
    )

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)


def test_get_parameter_dict():
    popt = [100, 10, 60, 0.6, 100, 0.2, 60]

    popt_dict = get_parameter_dict(popt)
    assert len(popt_dict) == 8
