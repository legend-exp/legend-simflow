from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from legendsimflow.drift_time import (
    drift_time_cost,
    drift_time_observables,
    remove_outliers,
)


@pytest.fixture
def rng():
    return np.random.default_rng(1)


def test_drift_time_observables(rng):
    x = rng.normal(1000, 100, 200_000)

    peak, q = drift_time_observables(x)
    assert peak == pytest.approx(1000, abs=5)
    assert q == pytest.approx(1000 + 100 * stats.norm.ppf(0.9), abs=2)
    assert drift_time_observables(x, percentile=50)[1] == pytest.approx(1000, abs=2)

    # NaNs are ignored
    x_nan = np.append(x, [np.nan] * 100)
    assert drift_time_observables(x_nan) == pytest.approx(drift_time_observables(x))

    # peak at the edge of the data
    assert drift_time_observables(np.full(100, 500.0)) == pytest.approx(
        [500] * 2, abs=1
    )


def test_drift_time_observables_peak_threshold(rng):
    # first peak at 2/3 of the height of the second
    x = np.concatenate([rng.normal(800, 50, 40_000), rng.normal(2000, 50, 60_000)])

    assert drift_time_observables(x)[0] == pytest.approx(800, abs=5)
    assert drift_time_observables(x, peak_threshold=0.8)[0] == pytest.approx(
        2000, abs=5
    )
    assert drift_time_observables(x, peak_threshold=1)[0] == pytest.approx(2000, abs=5)


def test_drift_time_observables_weights(rng):
    a = rng.normal(800, 50, 50_000)
    x = np.concatenate([a, rng.normal(2000, 50, 50_000)])
    w = np.concatenate([np.ones_like(a), np.zeros_like(a)])

    # zero weights remove the second component
    assert drift_time_observables(x, weights=w) == pytest.approx(
        drift_time_observables(a), abs=1
    )


def test_drift_time_cost(rng):
    x = rng.gamma(4, 150, 20_000)

    assert drift_time_cost(x, x) == 0
    # both observables move by the shift
    assert drift_time_cost(x, x + 10) == pytest.approx(2 * 10**2, rel=1e-2)


def test_remove_outliers():
    x = np.append(np.arange(1000.0), np.nan)

    assert remove_outliers(x, 99).max() == 989
    assert len(remove_outliers(x, 99)) == 990
    assert len(remove_outliers(x, 100)) == 1000
