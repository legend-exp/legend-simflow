from __future__ import annotations

import logging
import re
from pathlib import Path

import awkward as ak
import dbetto
import numpy as np
from dspeed.vis import WaveformBrowser
from lgdo import lh5
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import NDArray
from reboost.hpge.psd import _current_pulse_model as current_pulse_model
from scipy.optimize import curve_fit

from . import SimflowConfig, utils
from . import metadata as mutils

log = logging.getLogger(__name__)


def lookup_currmod_fit_data(
    hit_files: list[str, Path],
    lh5_group: str,
    ewin_center: float = 1593,
    ewin_width: float = 10,
) -> tuple[int, int]:
    """Extract the index of the event to fit.

    Considers events with |A/E| < 1.5 and finds the one that is closest to the
    median of the distribution.  Returns the index of the event in the file and
    the index of the file in the input file list.

    Parameters
    ----------
    hit_files
        tier-hit files used to determine the best index.
    lh5_group
        where the tier-hit data is found in the files.
    ewin_center
        center of the energy window to use for the event search (same units as
        in data).
    ewin_width
        width of the energy window to use for the event search (dame units as
        in data).
    """
    idxs = []
    energies = []
    dts = []

    for file in hit_files:
        energy = lh5.read(f"{lh5_group}/cuspEmax_ctc_cal", file).view_as("np")
        AoE = lh5.read(f"{lh5_group}/AoE_Classifier", file).view_as("np")

        # get drift time if possible
        if f"{lh5_group}/dt_eff" in lh5.ls(file):
            dt_eff = lh5.read(f"{lh5_group}/dt_eff", file).view_as("np")
        else:
            dt_eff = np.zeros(len(AoE))

        idx = np.where((abs(energy - ewin_center) < ewin_width / 2) & (abs(AoE) < 1.5))[
            0
        ]
        idxs.append(idx)
        energies.append(energy[idx])
        dts.append(dt_eff[idx])

        n = sum(len(d) for d in dts)

        if n > 500:
            break

    # now chose the best index (closest to median)
    if len(dts) == 0:
        msg = "no data found in the considered hit files"
        raise ValueError(msg)

    dts = ak.Array(dts)
    idxs = ak.Array(idxs)
    energies = ak.Array(energies)

    med = np.median(ak.flatten(dts))
    best = ak.argmin(abs(dts - med))
    n = ak.num(dts, axis=-1)

    array = np.full(len(ak.flatten(dts)), False)
    array[best] = True

    sel = ak.unflatten(array, n)

    idx = ak.flatten(idxs[sel])[0]
    file_idx = ak.where(ak.num(dts[sel], axis=-1) == 1)[0][0]

    return idx, file_idx


def fit_currmod(times: NDArray, current: NDArray) -> tuple:
    """Fit the model to the raw HPGe current pulse.

    Uses :func:`scipy.curve_fit` to fit
    :func:`reboost.hpge.psd._current_pulse_model` to the input raw pulse.

    Parameters
    ----------
    times
        the timesteps of the pulse.
    current
        the values of the pulse.

    Returns
    -------
        tuple of the best fit parameters, and arrays of the best fit model
        (time and current).
    """
    t = times
    A = current

    maxi = np.argmax(A)
    max_t = t[maxi]

    low = max_t - 1000
    high = max_t + 1000

    height = A[maxi]

    # initial guesses and ranges
    p0 = [height, max_t, 60, 0.6, 100, 0.2, 60]

    ranges = [
        [0, max_t - 100, 1, 0, 1, 0, 1],
        [height * 2, max_t + 100, 500, 1, 800, 1, 800],
    ]
    x = np.linspace(low, high, 1000)
    y = current_pulse_model(x, *p0)

    # do the fit
    popt, _ = curve_fit(
        current_pulse_model,
        t[(t > low) & (t < high)],
        A[(t > low) & (t < high)],
        p0=p0,
        bounds=ranges,
    )

    x = np.linspace(low, high, int(100 * (high - low)))
    y = current_pulse_model(x, *popt)

    return popt, x, y


def get_current_pulse(
    raw_file: Path | str,
    lh5_group: str,
    idx: int,
    dsp_config: str,
    dsp_output: str = "curr_av",
    align: str = "tp_aoe_max",
) -> tuple:
    """Extract the current pulse.

    Parameters
    ----------
    raw_file
        path to the raw tier file.
    lh5_group
        where to find the waveform table.
    idx
        the index of the waveform to read.
    dsp_config
        the :mod:`dspeed` configuration file defining the DSP processing chain
        to estimate the current pulse.
    dsp_output
        the name of the DSP output corresponding to the current pulse.
    align
        DSP value around which the pulses are aligned.
    """

    browser = WaveformBrowser(
        str(raw_file),
        lh5_group,
        dsp_config=dsp_config,
        lines=[dsp_output],
        align=align,
    )

    browser.find_entry(idx)

    t = browser.lines[dsp_output][0].get_xdata()
    A = browser.lines[dsp_output][0].get_ydata()

    return t, A


def plot_currmod_fit_result(
    t: NDArray, A: NDArray, model_t: NDArray, model_A: NDArray
) -> tuple[Figure, Axes]:
    """Plot the best fit results."""
    fig, ax = plt.subplots()

    ax.plot(t, A, linewidth=2, label="Waveform")
    ax.plot(model_t, model_A, color="red", linewidth=2, label="Fit")

    ax.legend()
    ax.set_xlim(-1200, 1200)

    ax.set_xlabel("Time [ns]")
    ax.set_ylabel("Current [ADC]")

    return fig, ax


def lookup_currmod_fit_inputs(
    config: SimflowConfig,
    runid: str,
    hpge: str,
    hit_tier_name: str = "hit",
) -> tuple[Path, int, Path]:
    """Find the raw file, event index and DSP configuration file.

    Parameters
    ----------
    config
        simflow configuration object.
    runid
        LEGEND-200 run identifier.
    hpge
        name of the HPGe detector
    hit_tier_name
        name of the hit tier. This is typically "hit" or "pht".
    """
    l200data = config.paths.l200data

    # get the reference cal run
    cal_runid = mutils.reference_cal_run(config.metadata, runid)
    _, period, run, _ = re.split(r"\W+", cal_runid)

    msg = f"inferred reference calibration run: {cal_runid}"
    log.debug(msg)

    # look for the dataflow config file
    df_cfgs = [p for p in l200data.glob("*config.*") if not p.name.startswith(".")]

    if len(df_cfgs) > 1:
        msg = f"found multiple configuration files in {l200data}, this cannot be!"
        raise RuntimeError(msg)

    msg = f"found dataflow configuration file: {df_cfgs[0]}"
    log.debug(msg)
    central_config = dbetto.utils.load_dict(df_cfgs[0])

    # get the paths to hit and raw tier files
    df_cfg = (
        central_config["setups"]["l200"]["paths"]
        if ("setups" in central_config)
        else central_config["paths"]
    )

    hit_path = Path(df_cfg[f"tier_{hit_tier_name}"].replace("$_", str(l200data)))
    msg = f"looking for hit tier files in {hit_path / 'cal' / period / run}/*"
    log.debug(msg)

    hit_files = list((hit_path / "cal" / period / run).glob("*"))

    if len(hit_files) == 0:
        msg = f"no hit tier files found in {hit_path}/cal/{period}/{run}"
        raise ValueError(msg)

    lh5_group = utils._get_lh5_table(config.metadata, hit_files[0], hpge, "hit", runid)

    msg = "looking for best event to fit"
    log.debug(msg)
    wf_idx, file_idx = lookup_currmod_fit_data(hit_files, lh5_group)

    raw_path = Path(df_cfg["tier_raw"].replace("$_", str(l200data)))
    hit_file = hit_files[file_idx]
    raw_file = raw_path / str(hit_file.relative_to(hit_path)).replace(
        hit_tier_name, "raw"
    )

    msg = f"determined raw file: {raw_file} (event index {wf_idx})"
    log.debug(msg)

    dsp_cfg_regex = r"l200-*-r%-T%-ICPC-dsp_proc_chain.*"
    dsp_cfg_files = list(
        (l200data / "inputs/dataprod/config/tier_dsp").glob(dsp_cfg_regex)
    )

    if dsp_cfg_files == []:
        dsp_cfg_files = list(
            (l200data / "inputs/dataprod/config/tier/dsp").glob(dsp_cfg_regex)
        )

    if len(dsp_cfg_files) != 1:
        msg = f"could not find a suitable dsp config file in {l200data} (or multiple found)"
        raise RuntimeError(msg)

    return raw_file.resolve(), wf_idx, dsp_cfg_files[0]
