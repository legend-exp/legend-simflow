# Copyright (C) 2026 Luigi Pertoldi <gipert@pm.me>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""Compare the simulated and data HPGe drift-time distributions.

The comparison uses the position of the first peak of the distribution and one
of its high percentiles.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks


def drift_time_observables(
    drift_times: ArrayLike,
    *,
    weights: ArrayLike | None = None,
    percentile: float = 90,
    smoothing: float = 50,
    peak_threshold: float = 0.5,
) -> NDArray:
    """First peak and a high percentile of a drift-time distribution.

    Both are read from the histogram of the drift times in 1 ns bins. The first
    peak is the lowest drift-time peak higher than `peak_threshold` times the
    global maximum of the histogram smoothed with a Gaussian of standard
    deviation `smoothing`. The percentile is the lower edge of the bin where the
    cumulative histogram reaches it.

    Parameters
    ----------
    drift_times
        drift time of each event in ns. NaNs are ignored.
    weights
        weight of each event, all equal by default.
    percentile
        percentile of the distribution, in percent.
    smoothing
        standard deviation of the Gaussian in ns.
    peak_threshold
        minimum peak height, as a fraction of the global maximum. ``1`` selects
        the global maximum.

    Returns
    -------
    obs
        ``[peak, percentile]`` in ns.
    """
    x = np.asarray(drift_times, dtype=float)
    finite = np.isfinite(x)
    x = x[finite]
    w = None if weights is None else np.asarray(weights, dtype=float)[finite]

    # padded so that a peak at the edge of the data is found too
    edges = np.arange(x.min() - smoothing, x.max() + smoothing + 1)
    hist = np.histogram(x, edges, weights=w)[0].astype(float)
    density = gaussian_filter1d(hist, smoothing, mode="constant")
    peak = edges[find_peaks(density, height=peak_threshold * density.max())[0][0]]

    cdf = np.cumsum(hist)
    q = edges[np.searchsorted(cdf, percentile / 100 * cdf[-1])]
    return np.array([peak + 0.5, q])


def drift_time_cost(
    data_drift_times: ArrayLike,
    sim_drift_times: ArrayLike,
    sim_weights: ArrayLike | None = None,
    **obs_kwargs,
) -> float:
    """Distance between the data and simulated drift-time distributions.

    Sum of the squared differences (ns²) between the
    :func:`drift_time_observables` of the simulation and of the data.

    Parameters
    ----------
    data_drift_times
        drift time of each data event in ns.
    sim_drift_times
        drift time of each simulated event in ns.
    sim_weights
        weight of each simulated event, all equal by default.
    **obs_kwargs
        passed to :func:`drift_time_observables`.
    """
    diff = drift_time_observables(sim_drift_times, weights=sim_weights, **obs_kwargs)
    diff -= drift_time_observables(data_drift_times, **obs_kwargs)
    return float(np.sum(diff**2))


def remove_outliers(values: ArrayLike, percentile: float = 99) -> NDArray:
    """Values up to the `percentile` (in percent) of the sample; NaNs are dropped."""
    x = np.asarray(values)
    return x[x <= np.nanpercentile(x, percentile)]
