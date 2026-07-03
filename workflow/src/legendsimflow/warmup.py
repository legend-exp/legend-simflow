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
"""Serial warmup of heavy dependencies and Numba caches.

Run once before the parallel Snakemake jobs (via ``python -m
legendsimflow.warmup``, wired to the ``warmup`` pixi task). A lazily-compiled
``@njit(cache=True)`` kernel is compiled and written to a shared on-disk cache
only on its first *call*; parallel jobs racing to write it segfault. Warming
each kernel once here leaves them only reading the cache.
"""
# heavy imports are for their side effects, all inside the warmup function
# ruff: noqa: F401, ICN001, PLC0415

from __future__ import annotations


def warm_numba_caches() -> None:
    """Warm the heavy imports and the Numba kernels the workflow calls."""
    # importing these pays their one-off cost (Matplotlib font cache, ...) and
    # compiles their vectorized kernels (dspeed/lh5 were ``import *``, illegal
    # inside a function; a plain import triggers the same compilation)
    import dspeed.processors
    import lh5.compression
    import matplotlib
    import numpy as np
    import pygama
    import pygeomhpges
    import pygeomoptics
    import pygeomtools
    import reboost
    import remage
    import revertex

    # lazily-compiled kernels must be called once, with the exact runtime type
    # signature (Numba caches per signature). These are the only cache=True
    # kernels the parallel rules reach.
    from reboost.hpge.psd import _current_pulse_model

    from .reboost import _cluster_photoelectrons_flat

    # currmod fit: float64 time array + seven float64 scalars
    _current_pulse_model(
        np.linspace(-500.0, 1500.0, 16), 1.0, 0.0, 60.0, 0.6, 100.0, 0.2, 60.0
    )
    # opt tier: int64 list offsets, float64 time/amplitude content
    _cluster_photoelectrons_flat(
        np.array([0, 2, 3], dtype=np.int64),
        np.array([0.0, 1.0, 5.0]),
        np.array([1.0, 2.0, 3.0]),
        10.0,
    )


if __name__ == "__main__":
    warm_numba_caches()
