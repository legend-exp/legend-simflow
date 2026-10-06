from __future__ import annotations

import pytest

from legendsimflow.impurity_tuning import get_wf_chi2


def test_get_wf_chi2():
    elecmod = {
        "slope_0": {"dep_0": {"rms": 0.1}, "dep_2": {"rms": 0.2}},
        "slope_3": {"dep_1": {"rms": 0.4}},
    }
    grid_info = {"slope_min": -1, "slope_step": 0.5, "dep_min": 1000, "dep_step": 100}

    dep, slope, chi2 = get_wf_chi2(elecmod, grid_info, wf_scale=0.1)

    assert dep == [1000, 1200, 1100]
    assert slope == [-1, -1, 0.5]
    assert chi2 == pytest.approx([1, 4, 16])
