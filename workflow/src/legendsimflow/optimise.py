# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>
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

"""Utilities to compute optimal simulation settings from benchmark data.

Given a measured simulation speed (events/second) and user-defined
constraints, the functions in this module compute recommended values for
``primaries_per_job`` and ``number_of_jobs`` to be placed in the Simflow
configuration files.

The algorithm:

1. **Maximise** ``primaries_per_job`` within the ``max_runtime_mins``
   constraint, then round *down* to ``n_sig_figs`` significant figures so
   that jobs are guaranteed to stay under the time limit.
2. **Minimise** ``number_of_jobs`` so that the product
   ``primaries_per_job * number_of_jobs >= target_total_primaries``, then
   round *up* to ``n_sig_figs`` significant figures so the statistics
   requirement is guaranteed to be met.
"""

from __future__ import annotations

import math


def round_down_nsf(x: float, n: int = 2) -> int:
    """Round *down* to ``n`` significant figures.

    For values whose order of magnitude is already smaller than ``n``
    significant figures can express (i.e. ``x < 10^(n-1)``), the value is
    simply truncated to an integer.

    Parameters
    ----------
    x
        Input value (must be positive).
    n
        Number of significant figures to keep.

    Returns
    -------
    int
        The value rounded down to ``n`` significant figures.

    Examples
    --------
    >>> round_down_nsf(1291860, 2)
    1200000
    >>> round_down_nsf(717.7, 2)
    710
    >>> round_down_nsf(3.6, 2)
    3
    """
    if x <= 0:
        return 0
    k = math.floor(math.log10(x))
    magnitude = k - n + 1
    if magnitude <= 0:
        # x already has at most n significant figures; just truncate to int
        return max(1, int(x))
    factor = 10**magnitude
    return math.floor(x / factor) * factor


def round_up_nsf(x: float, n: int = 2) -> int:
    """Round *up* to ``n`` significant figures.

    For values whose order of magnitude is already smaller than ``n``
    significant figures can express (i.e. ``x < 10^(n-1)``), the value is
    simply rounded up to the next integer.

    Parameters
    ----------
    x
        Input value (must be positive).
    n
        Number of significant figures to keep.

    Returns
    -------
    int
        The value rounded up to ``n`` significant figures.

    Examples
    --------
    >>> round_up_nsf(314, 2)
    320
    >>> round_up_nsf(31, 2)
    31
    >>> round_up_nsf(3.7, 2)
    4
    """
    if x <= 0:
        return 0
    k = math.floor(math.log10(x))
    magnitude = k - n + 1
    if magnitude <= 0:
        # x already has at most n significant figures; just take the ceiling
        return max(1, math.ceil(x))
    factor = 10**magnitude
    return math.ceil(x / factor) * factor


def compute_primaries_per_job(
    speed: float,
    max_runtime_mins: float,
    n_sig_figs: int = 2,
) -> int:
    """Compute the number of primaries per job satisfying the runtime constraint.

    Computes the maximum number of events that can be simulated within
    ``max_runtime_mins`` at the given ``speed``, then rounds *down* to
    ``n_sig_figs`` significant figures.

    Parameters
    ----------
    speed
        Simulation speed in events/second, as measured from a benchmark run.
    max_runtime_mins
        Maximum allowed job runtime in minutes.
    n_sig_figs
        Number of significant figures for rounding (default: 2).

    Returns
    -------
    int
        Recommended ``primaries_per_job`` value.
    """
    raw = speed * max_runtime_mins * 60
    return max(1, round_down_nsf(raw, n_sig_figs))


def compute_number_of_jobs(
    primaries_per_job: int,
    target_total_primaries: int,
    n_sig_figs: int = 2,
) -> int:
    """Compute the minimum number of jobs to meet the total primaries target.

    Divides ``target_total_primaries`` by ``primaries_per_job``, then rounds
    *up* to ``n_sig_figs`` significant figures to guarantee that the total
    simulated statistics requirement is met.

    Parameters
    ----------
    primaries_per_job
        Number of primaries simulated in each individual job.
    target_total_primaries
        Minimum number of total primary events required.
    n_sig_figs
        Number of significant figures for rounding (default: 2).

    Returns
    -------
    int
        Recommended ``number_of_jobs`` value.
    """
    raw = math.ceil(target_total_primaries / primaries_per_job)
    return max(1, round_up_nsf(raw, n_sig_figs))


def compute_sim_settings(
    speed: float,
    max_runtime_mins: float,
    target_total_primaries: int,
    n_sig_figs: int = 2,
) -> tuple[int, int]:
    """Compute optimal ``primaries_per_job`` and ``number_of_jobs``.

    Given the simulation speed from a benchmark run and the user-defined
    constraints, returns the recommended values for ``primaries_per_job`` and
    ``number_of_jobs`` to insert into the Simflow configuration.

    The following constraints are satisfied:

    * Each job runs for at most ``max_runtime_mins`` minutes.
    * The total simulated statistics exceed ``target_total_primaries``.
    * The total number of jobs is minimised.
    * Both values have at most ``n_sig_figs`` significant figures.

    Parameters
    ----------
    speed
        Simulation speed in events/second, as measured from a benchmark run.
    max_runtime_mins
        Maximum allowed job runtime in minutes.
    target_total_primaries
        Minimum number of total primary events required.
    n_sig_figs
        Number of significant figures for rounding (default: 2).

    Returns
    -------
    tuple[int, int]
        A tuple ``(primaries_per_job, number_of_jobs)``.

    Examples
    --------
    >>> compute_sim_settings(speed=500, max_runtime_mins=30, target_total_primaries=1_000_000)
    (900000, 2)
    """
    ppj = compute_primaries_per_job(speed, max_runtime_mins, n_sig_figs)
    njobs = compute_number_of_jobs(ppj, target_total_primaries, n_sig_figs)
    return ppj, njobs
