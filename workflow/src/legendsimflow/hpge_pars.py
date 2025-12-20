from __future__ import annotations

import inspect
from pathlib import Path

import awkward as ak
import dbetto
import numpy as np
from dspeed.vis import WaveformBrowser
from legendmeta import LegendMetadata
from lgdo import lh5
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import ArrayLike, NDArray
from reboost.hpge.psd import _current_pulse_model as model
from scipy.optimize import curve_fit


def get_index(
    files: list, lh5_group: str, peak: float = 1593, peak_range: float = 5
) -> tuple[int, int]:
    """Extract the index of the event to fit"""
    idxs = []
    energies = []
    dts = []

    for file in files:
        energy = lh5.read(f"{lh5_group}/cuspEmax_ctc_cal", file).view_as("np")
        AoE = lh5.read(f"{lh5_group}/AoE_Classifier", file).view_as("np")

        # get drift time if possible
        if f"{lh5_group}/dt_eff" in lh5.ls(file):
            dt_eff = lh5.read(f"{lh5_group}/dt_eff", file).view_as("np")
        else:
            dt_eff = np.zeros(len(AoE))

        idx = np.where((abs(energy - peak) < peak_range) & (abs(AoE) < 1.5))[0]
        idxs.append(idx)
        energies.append(energy[idx])
        dts.append(dt_eff[idx])
        n = sum(len(d) for d in dts)

        if n > 500:
            break

    # now chose the best index (closest to median)

    if len(dts) == 0:
        msg = "No data found in files"
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


def get_raw_file(
    hit_file: Path, hit_path: Path, raw_path: Path, hit_name: str = "pht"
) -> Path:
    """Get the raw file from the paths"""

    raw_file = hit_file.relative_to(hit_path)
    file = Path(raw_path) / raw_file

    return file.parent / file.name.replace(hit_name, "raw")


def fit_waveform(times: NDArray, current: NDArray) -> tuple:
    """Fit the waveform

    Parameters
    ----------
    times
        The timesteps of the waveform
    current
        The values of the waveform

    Returns
    -------
        tuple of the best fit parameters, and arrays of the best fit model (time and current).
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
    y = model(x, *p0)

    # do the fit
    popt, _ = curve_fit(
        model,
        t[(t > low) & (t < high)],
        A[(t > low) & (t < high)],
        p0=p0,
        bounds=ranges,
    )

    x = np.linspace(low, high, int(100 * (high - low)))
    y = model(x, *popt)

    return popt, x, y


def get_waveform(
    raw_file: Path | str,
    lh5_group: str,
    idx: int,
    dsp_config: str,
    *,
    current: str = "curr_av",
    align: str = "tp_aoe_max",
) -> tuple:
    """Extract the current waveform"""

    browser = WaveformBrowser(
        str(raw_file),
        lh5_group,
        dsp_config=dsp_config,
        lines=[current],
        align=align,
    )

    browser.find_entry(idx)

    t = browser.lines[current][0].get_xdata()
    A = browser.lines[current][0].get_ydata()

    return t, A


def plot_fit(
    t: NDArray, A: NDArray, model_t: NDArray, model_A: NDArray
) -> tuple[Figure, Axes]:
    """Plot the best fit reconstruction"""
    fig, ax = plt.subplots()

    ax.plot(t, A, linewidth=2, label="Waveform")
    ax.plot(model_t, model_A, color="red", linewidth=2, label="Fit")

    ax.legend()
    ax.set_xlim(-1200, 1200)

    ax.set_xlabel("time [ns]")
    ax.set_ylabel("Current [ADC]")

    return fig, ax


def get_parameter_dict(popt: ArrayLike) -> dict:
    """Get the parameter results as a dictionary"""

    params = list(inspect.signature(model).parameters)
    param_names = params[1:]

    popt_dict = dict(zip(param_names, popt, strict=True))
    popt_dict["mean_AoE"] = popt_dict["amax"] / 1593

    for key, value in popt_dict.items():
        popt_dict[key] = float(f"{value:.3g}")

    return popt_dict


def get_pulse_pars(
    hpge: str,
    prod_cycle_path: str,
    period: str,
    run: str,
    plot_path: str | None = None,
    hit_tier: str = "pht",
    use_channel_name: bool = False,
    meta: LegendMetadata | None = None,
) -> dict:
    """Get the parameters of the fitted waveform

    Parameters
    ----------
    hpge
        The hpgeector name
    prod_cycle_path
        The path to the production cycle with data.
    period
        The period to inspect.
    run
        The run to inspect.
    plot_path
        Optional path to save a plot.
    hit_tier
        The name of the hit level tier (typically hit or pht)
    use_channel_name
        Whether the data is the new format referenced by channel name, or the old one.
    meta
        The legendMetadata instance, if `None` will be constructed.
    """

    # extract some metadata

    config_path = next(
        p
        for p in Path(f"{prod_cycle_path}/").glob("*config.*")
        if not p.name.startswith(".")
    )
    central_config = dbetto.utils.load_dict(str(config_path))

    if meta is None:
        meta = LegendMetadata()

    timestamp = meta.datasets.runinfo[period][run].cal.start_key

    chmap = meta.channelmap(timestamp)
    rawid = chmap[hpge].daq.rawid

    # get the paths
    c_dict = (
        central_config["setups"]["l200"]["paths"]
        if ("setups" in central_config)
        else central_config["paths"]
    )

    hit_path = c_dict[f"tier_{hit_tier}"].replace("$_", prod_cycle_path)

    raw_path = c_dict["tier_raw"].replace("$_", prod_cycle_path)

    files = list(Path(f"{hit_path}/cal/{period}/{run}/").glob("*"))

    if len(files) == 0:
        msg = f"No files found in {hit_path}/cal/{period}/{run}"
        raise ValueError(msg)

    idx, file_idx = get_index(
        files, f"hit/{hpge}" if use_channel_name else f"ch{rawid}/hit"
    )

    raw_file = get_raw_file(
        files[file_idx], hit_path, raw_path, hit_name=hit_tier
    ).resolve()

    # extract the waveform to fit

    if Path(f"{prod_cycle_path}/inputs/dataprod/config/tier_dsp/").exists():
        path = Path(f"{prod_cycle_path}/inputs/dataprod/config/tier_dsp/")
    elif Path(f"{prod_cycle_path}/inputs/dataprod/config/tier/dsp/").exists():
        path = Path(f"{prod_cycle_path}/inputs/dataprod/config/tier/dsp/")
    else:
        msg = "No path with dsp config file found!"
        raise RuntimeError(msg)

    config = str(next(path.glob("l200-*-r%-T%-ICPC-dsp_proc_chain.*")))

    t, A = get_waveform(
        raw_file,
        current="curr_av",
        lh5_group=f"raw/{hpge}" if use_channel_name else f"ch{rawid}/raw",
        dsp_config=config,
        idx=idx,
    )

    # now perform the fit
    popt, x, y = fit_waveform(t, A)

    # now plot
    if plot_path is not None:
        fig, _ = plot_fit(t, A, x, y)
        plt.savefig(f"{plot_path}/current_fit-{period}-{run}-{hpge}.pdf")
        plt.close(fig)

    # get the results as a dict
    popt_dict = get_parameter_dict(popt)

    return popt_dict, t, A
