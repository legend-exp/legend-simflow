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
legendsimflow.warmup [simflow-config.yaml]``, wired to the ``warmup`` pixi
task). A lazily-compiled ``@njit(cache=True)`` kernel is compiled and written to
a shared on-disk cache only on its first *call*; parallel jobs racing to write
it segfault. Warming each kernel once here leaves them only reading the cache.
"""
# heavy imports are kept inside the function so importing this module stays
# cheap (e.g. sphinx autodoc can import it without the full runtime stack).
# ruff: noqa: F401, ICN001, PLC0415

from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger(__name__)


def warm_numba_caches(simflow_config: str | Path | dict | None = None) -> None:
    """Warm the heavy imports and the Numba kernels the workflow calls."""
    # importing these pays their one-off cost (Matplotlib font cache, ...) and
    # compiles their vectorized kernels
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

    # dspeed.processors/lh5.compression expose their members lazily (PEP 562), so
    # a plain import leaves their eager @guvectorize kernels (e.g.
    # moving_window_multi, used by the parallel realistic-PSL rule) uncompiled and
    # the parallel jobs race the cache and segfault. The precompile task this
    # replaced used `import *` to force every name to load and compile; `import *`
    # is illegal inside a function, so reproduce it by materializing __all__.
    for _module in (dspeed.processors, lh5.compression):
        for _name in _module.__all__:
            getattr(_module, _name)

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

    # HPGe DSP-chain kernels raced by the parallel par-tier rules, warmed only
    # when this production actually runs them (see _resolve_dsp_config). Warn and
    # skip on any problem rather than aborting: warmup gates every rule, and the
    # par-tier rules would surface a real misconfiguration themselves.
    try:
        dsp_config = _resolve_dsp_config(simflow_config)
        if dsp_config is not None:
            log.info("warming the HPGe DSP-chain with config %s", dsp_config)
            warm_hpge_dsp_cache(dsp_config)
    except Exception:
        log.warning(
            "could not warm the HPGe DSP-chain; the parallel par-tier rules may "
            "race the Numba cache and segfault",
            exc_info=True,
        )


def _resolve_dsp_config(simflow_config: str | Path | dict | None) -> Path | None:
    """DSP config to warm, or ``None`` if the par-tier rules that race it won't run."""
    if simflow_config is None:
        return None
    if isinstance(simflow_config, (str, Path)) and not Path(simflow_config).exists():
        return None

    import dbetto

    from legendsimflow import utils
    from legendsimflow.metadata import get_tier_settings

    # cheap peek before the heavier full context init
    raw = (
        simflow_config
        if isinstance(simflow_config, dict)
        else dbetto.utils.load_dict(str(simflow_config))
    )
    if "l200data" not in raw.get("paths", {}):
        return None

    config = utils.init_simflow_context(simflow_config).config
    hit = get_tier_settings(config, "hit")
    if not (hit.get("simulate_psd", True) and hit.get("simulate_psd_with_psl", True)):
        return None

    return utils.lookup_dsp_config(config.paths.l200data)


def _make_dummy_raw_file(path: str, lh5_group: str = "ch0/raw", n_wfs: int = 4) -> str:
    """Write a small synthetic HPGe raw file to drive the DSP-chain warmup, return the LH5 group written."""
    import lh5
    import numpy as np
    from lgdo import Array, Table, WaveformTable
    from reboost.hpge.psd import _current_pulse_model

    # shaped current -> integrated charge (see tests/conftest.py make_cal_data)
    t = np.linspace(0, 99968, 6249)
    curr = _current_pulse_model(t, 1.0, 51000.0, 50.0, 0.55, 150.0, 0.2, 80.0)
    curr /= np.sum(curr) * 16
    charge = np.cumsum(curr)
    windowed = charge[2625:4025]  # 1400-sample high-resolution window
    presummed = charge[::8][:-1] * 8  # 781-sample presummed (16 ns -> 128 ns)

    energy = np.full(n_wfs, 1593.0)
    win = (np.vstack([windowed] * n_wfs) * energy[:, None]).astype("int32")
    presum = (np.vstack([presummed] * n_wfs) * energy[:, None]).astype("int32")

    tb = Table(size=n_wfs)
    tb.add_field(
        "waveform_presummed",
        WaveformTable(values=presum, dt=128, dt_units="ns", t0=0, t0_units="ns"),
    )
    tb.add_field(
        "waveform_windowed",
        WaveformTable(values=win, dt=16, dt_units="ns", t0=42000, t0_units="ns"),
    )
    tb.add_field("presum_rate", Array(np.full(n_wfs, 8, dtype="float32")))

    lh5.write(tb, lh5_group, path, wo_mode="of")
    return lh5_group


def warm_hpge_dsp_cache(dsp_config: str | Path) -> None:
    """Compile the HPGe DSP-chain (and fit) kernels raced by the par-tier rules."""
    import tempfile
    import warnings

    import numpy as np
    from scipy.stats import norm

    from legendsimflow import hpge_pars, superpulses

    dsp = str(dsp_config)
    # only compiling kernels here; numpy warnings from the synthetic pulse are noise
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        with tempfile.TemporaryDirectory() as tmp:
            raw_file = str(Path(tmp) / "dummy-raw.lh5")
            lh5_group = _make_dummy_raw_file(raw_file)
            n = 4

            # current-pulse chain (also warms the tp_aoe_max alignment conversions)
            hpge_pars.get_current_pulse(
                raw_file, lh5_group, 0, dsp, "curr_av", "tp_aoe_max"
            )
            # superpulse chain: charge (wf_pz_win) + current + baseline + energy outputs
            superpulses.get_wfs_for_slice(
                [raw_file], lh5_group, list(range(n)), [0] * n, dsp_config=dsp
            )

        # iminuit noise fit (deterministic Gaussian sample, no RNG)
        hpge_pars.fit_noise_gauss(norm.ppf(np.linspace(1e-3, 1 - 1e-3, 2000)), bins=100)


if __name__ == "__main__":
    import sys

    # the prod invocation runs from the prod-cycle cwd, where the config lives;
    # override the path as first arg (e.g. for the test tasks)
    config_file = sys.argv[1] if len(sys.argv) > 1 else "simflow-config.yaml"
    warm_numba_caches(config_file)
