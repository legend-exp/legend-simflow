from __future__ import annotations

import functools
import logging
import re
from collections.abc import Callable
from pathlib import Path
from typing import Any

import awkward as ak
import numpy as np
from dbetto import AttrsDict, TextDB
from iminuit import Minuit, cost
from legendmeta import LegendMetadata
from lgdo import lh5
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import ArrayLike, NDArray
from pygama.math.distributions import gaussian
from reboost.hpge.psd import _current_pulse_model as current_pulse_model
from scipy.optimize import curve_fit

from . import metadata as mutils
from . import utils

log = logging.getLogger(__name__)


def lookup_currmod_fit_data(
    hit_files: list[str, Path],
    lh5_group: str,
    ewin_center: float = 1593,
    ewin_width: float = 10,
) -> tuple[int, int]:
    """Extract the index of the event to fit.

    Considers events with ``abs(A/E) < 1.5`` and finds the one that is closest
    to the median of the distribution.  Returns the index of the event in the
    file and the index of the file in the input file list.

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
        aoe = lh5.read(f"{lh5_group}/AoE_Classifier", file).view_as("np")

        # get drift time if possible
        if f"{lh5_group}/dt_eff" in lh5.ls(file):
            dt_eff = lh5.read(f"{lh5_group}/dt_eff", file).view_as("np")
        else:
            dt_eff = np.zeros(len(aoe))

        idx = np.where((abs(energy - ewin_center) < ewin_width / 2) & (abs(aoe) < 1.5))[
            0
        ]
        idxs.append(idx)
        energies.append(energy[idx])
        dts.append(dt_eff[idx])

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


def fit_noise_gauss(
    data: ArrayLike,
    bins: int,
    *,
    fit_range: tuple | None = None,
    sigma_range: tuple | None = None,
) -> Minuit:
    """Fit the data to a Gaussian to extract the resolution.

    Performs a binned maximum likelihood fit using `minuit`.

    Parameters
    ----------
    data
        an array of the data to fit.
    bins
        The number of bins.
    fit_result
        The results of the `iminuit` fit.
    fit_range
        The range to use for the fit, if `None` this is determined from the data as +/- 5 standard deviations round the mean.
    sigma_range
        The range of sigma values for the fit, if `None` is determined from the data.

    Returns
    -------
        The minuit object holding the fit results.
    """

    if fit_range is None:
        fit_range = (np.mean(data) - 5 * np.std(data), np.mean(data) + 5 * np.std(data))

    if sigma_range is None:
        sigma_range = (0.1 * np.std(data), 10 * np.std(data))

    n, xe = np.histogram(data, bins=bins, range=fit_range)
    c = cost.BinnedNLL(n, xe, gaussian.cdf_norm)

    # initialise the fit
    m_norm = Minuit(
        c,
        x_lo=fit_range[0],
        x_hi=fit_range[1],
        mu=(fit_range[0] + fit_range[1]) / 2,
        sigma=(sigma_range[0] + sigma_range[1]) / 2,
    )

    # set parameters
    m_norm.fixed["x_lo", "x_hi"] = True
    m_norm.limits["mu"] = fit_range
    m_norm.limits["sigma"] = sigma_range

    # perform the fit
    m_norm.migrad()
    m_norm.hesse()

    return m_norm


def plot_noise_waveforms(
    noise: ArrayLike, temp: ArrayLike, norm: float = 1
) -> tuple[Figure, Any]:
    """Plot the waveforms with noise and the noise alone."""

    temp = norm * temp / np.max(temp)

    fig, axs = plt.subplots(2, 1, figsize=(6, 4), sharex=True)

    for idx, wf in enumerate(noise[0:20]):
        axs[0].plot(wf + temp, label=("noise + current template" if idx == 0 else None))
        axs[1].plot(wf, label=("noise" if idx == 0 else None))

    axs[1].set_xlabel("time [sample]")

    axs[0].set_ylabel("waveform")
    axs[1].set_ylabel("waveform")

    axs[0].legend()
    axs[1].legend()

    return fig, axs


def plot_gauss_fit(
    data: ArrayLike,
    fit_result: Minuit,
    fit_range: tuple | None = None,
    bins: int = 100,
    nominal_val: float | None = None,
) -> tuple[Figure, Axes]:
    """Plot the result of the Gaussian fit.

    Parameters
    ----------
    data
        an array of the data to fit.
    bins
        The number of bins.
    fit_range
        The range to use for the fit, if `None` this is determined from the data as +/- 5 standard deviations round the mean.
    nominal_val
        The nominal mean to add as a line on the plot.

    """
    if fit_range is None:
        fit_range = (np.mean(data) - 5 * np.std(data), np.mean(data) + 5 * np.std(data))

    fig, ax = plt.subplots(figsize=(6, 4))

    ax.hist(data, bins=bins, range=fit_range, alpha=0.8, density=True)

    x = np.linspace(*fit_range, 10000)

    ax.plot(
        x,
        gaussian.pdf_norm(
            x,
            x_lo=fit_result.values["x_lo"],
            x_hi=fit_result.values["x_hi"],
            mu=fit_result.values["mu"],
            sigma=fit_result.values["sigma"],
        ),
        label="Fit",
    )

    if nominal_val is not None:
        ax.axvline(nominal_val, color="red", linestyle="--", label="Nominal value")

    ax.set_xlabel("A$_{max}$ [arb]")
    ax.set_ylabel("Prob [arb.]")
    ax.legend()

    ax.set_yscale("linear")
    ax.get_legend().set_title(r"$\sigma$" + f" = {fit_result.values['sigma']:.2f}")

    return fig, ax


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
    # HACK: importing this messes up pint registries
    from dspeed.vis import WaveformBrowser  # noqa: PLC0415

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


def get_noise_waveforms(
    raw_files: list,
    hit_files: list,
    lh5_group: str,
    dsp_config: str,
    dsp_output: str,
    *,
    threshold: float = 5,
    length: int = 1000,
    maximum_number: int | None = None,
    energy_var: str = "cuspEmax_cal",
):
    """Extract a matrix of noise waveforms, after applying some DSP processing.

    The waveforms are only those with energy less than `threshold` and are the
    result of the dsp processing defined in `config_path` with output variable
    `dsp_output`.

    Parameters
    ----------
    raw_files
        List of paths to raw files.
    hit_files
        List of paths to hit files.
    lh5_group
        The name of the lh5_group to find the waveform table in.
    dsp_config
        the :mod:`dspeed` configuration file defining the DSP processing chain
        to estimate the current pulse.
    dsp_output
        the name of the DSP output corresponding to the current pulse.
    threshold
        energy threshold to apply to select the noise waveforms
    maximum_number
        Number of waveforms to extract.
    length
        length of noise waveform to extract (takes the first `N` samples).

    Returns
    -------
    a 2D array of the waveforms.

    """
    from dspeed.vis import WaveformBrowser  # noqa: PLC0415

    waveforms = []

    for raw_file, hit_file in zip(raw_files, hit_files, strict=True):
        browser = WaveformBrowser(
            str(raw_file), lh5_group, dsp_config=dsp_config, lines=[dsp_output]
        )
        energies = lh5.read(
            f"{lh5_group.replace('raw', 'hit')}/{energy_var}", hit_file
        ).view_as("np")

        # only look at low energy
        indices = np.where(energies < threshold)

        for idx in indices[0]:
            browser.find_entry(idx, append=False)

            waveform = browser.lines[dsp_output][0].get_ydata()
            waveforms.append(waveform[:length].reshape(1, length))

            if maximum_number is not None and len(waveforms) >= maximum_number:
                return np.concatenate(waveforms, axis=0)

    return np.concatenate(waveforms, axis=0)


def plot_currmod_fit_result(
    t: NDArray, A: NDArray, model_t: NDArray, model_A: NDArray
) -> tuple[Figure, Axes]:
    """Plot the best fit results."""
    fig, ax = plt.subplots(figsize=(6, 4))

    ax.plot(
        t,
        A,
        linewidth=2,
        label="Current signal",
    )
    ax.plot(model_t, model_A, label="Model", color="tab:red")

    ax.legend()

    ax.set_xlabel("Time [ns]")
    ax.set_ylabel("Current [ADC]")

    return fig, ax


def estimate_mean_aoe(popt: list, energy: float = 1593) -> float:
    """Estimate the maximum aoe from the parameters of the `current_pulse_model` `popt`."""
    # get the maximum of the template
    x = np.linspace(-1000, 3000, 4001)
    temp = current_pulse_model(x, *popt)

    return float(np.max(temp) / energy)


def get_waveform_maxima(
    template: ArrayLike, noise_wfs: ArrayLike, *, norm: float = 1
) -> NDArray:
    """Extract the maximum of each waveform based on combining the
    `template` with each waveform in `noise_wfs`.

    Note
    ----
    The length of the template must be the same as the waveforms in `noise_wfs`

    Parameters
    ----------
    template
        The template of the waveform to use.
    noise_wfs
        2D array of each noise waveform.
    norm
        The normalisation for the template.
    """
    # normalise the template
    template = norm * template / np.max(template)

    maximum = np.max(noise_wfs + template, axis=1)

    # remove nan
    maximum = maximum[~np.isnan(maximum)]

    # remove far outliers
    return maximum[abs(maximum - np.mean(maximum)) < np.std(maximum) * 5]


def lookup_file_paths(l200data: str, runid: str, hit_tier_name: str) -> AttrsDict:
    """Lookup the paths to the `hit` and `raw` files."""

    _, period, run, data_type = re.split(r"\W+", runid)

    if isinstance(l200data, str):
        l200data = Path(l200data)

    df_cfg = utils.lookup_dataflow_config(l200data).paths

    hit_path = df_cfg[f"tier_{hit_tier_name}"]
    hit_files = list((hit_path / data_type / period / run).glob("*.lh5"))

    raw_files = [
        Path(
            str(hit_file)
            .replace(str(hit_path), str(df_cfg.tier_raw))
            .replace(hit_tier_name, "raw")
        )
        for hit_file in hit_files
    ]

    return AttrsDict({"raw": raw_files, "hit": hit_files})


def lookup_currmod_fit_inputs(
    l200data: str | Path,
    metadata: LegendMetadata,
    runid: str,
    hpge: str,
    hit_tier_name: str = "hit",
) -> tuple[Path, int, Path]:
    """Find the raw file, event index and DSP configuration file.

    Parameters
    ----------
    l200data
        The path to the L200 data production cycle.
    metadata
        The metadata instance
    runid
        LEGEND-200 run identifier, must be of the form `{EXPERIMENT}-{PERIOD}-{RUN}-{TYPE}`.
    hpge
        name of the HPGe detector
    hit_tier_name
        name of the hit tier. This is typically "hit" or "pht".
    """
    if isinstance(l200data, str):
        l200data = Path(l200data)

    df_cfg = utils.lookup_dataflow_config(l200data).paths

    # get the reference cal run
    cal_runid = mutils.reference_cal_run(metadata, runid)
    _, period, run, _ = re.split(r"\W+", cal_runid)

    msg = f"inferred reference calibration run: {cal_runid}"
    log.debug(msg)

    hit_path = df_cfg[f"tier_{hit_tier_name}"]
    msg = f"looking for hit tier files in {hit_path / 'cal' / period / run}/*"
    log.debug(msg)

    hit_files = list((hit_path / "cal" / period / run).glob("*"))

    if len(hit_files) == 0:
        msg = f"no hit tier files found in {hit_path}/cal/{period}/{run}"
        raise ValueError(msg)

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

    lh5_group = mutils._get_lh5_table(metadata, hit_files[0], hpge, "hit", runid)

    msg = "looking for best event to fit"
    log.debug(msg)

    wf_idx, file_idx = lookup_currmod_fit_data(hit_files, lh5_group)

    raw_path = df_cfg.tier_raw
    hit_file = hit_files[file_idx]
    raw_file = (
        raw_path / str(hit_file.relative_to(hit_path)).replace(hit_tier_name, "raw")
    ).resolve()
    msg = f"determined raw file: {raw_file} (event index {wf_idx})"
    log.debug(msg)

    return raw_file, wf_idx, dsp_cfg_files[0]


def lookup_energy_res_metadata(
    l200data: str | Path,
    metadata: LegendMetadata,
    runid: str,
    *,
    hit_tier_name: str = "hit",
    pars_db: TextDB | None = None,
) -> AttrsDict:
    r"""Lookup the measured HPGe energy resolution metadata from LEGEND-200 data.

    The metadata refers to the following model:

    .. math::

        \text{FWHM}(E) = \sqrt{a + bE}

    where :math:`E` is in keV.

    Returns
    -------
    Mapping of HPGe name to metadata dictionary.

    Parameters
    ----------
    l200data
        The path to the L200 data production cycle.
    metadata
        The metadata instance
    runid
        LEGEND-200 run identifier, must be of the form `{EXPERIMENT}-{PERIOD}-{RUN}-{TYPE}`.
    hit_tier_name
        name of the hit tier. This is typically "hit" or "pht".
    pars_db
        optional existing *non-lazy* instance of
        ``TextDB(".../path/to/prod/generated/par_{hit_tier_name}")``.
    """
    if hit_tier_name not in ("hit", "pht"):
        raise NotImplementedError

    if isinstance(l200data, str):
        l200data = Path(l200data)

    # get the paths to generated parameters
    if pars_db is None:
        pars_db = utils.init_generated_pars_db(l200data, tier=hit_tier_name, lazy=True)

    msg = f"loading {hit_tier_name} pars of production {l200data}"
    log.debug(msg)

    # get the pars file at the correct timestamp
    tstamp = mutils.runinfo(metadata, runid).start_key
    chmap = metadata.hardware.configuration.channelmaps.on(tstamp)
    pars_file = pars_db.on(tstamp)

    out_dict = {}
    for key, detmeta in pars_file.items():
        # handle data prod formats
        hpge = (
            chmap.map("daq.rawid")[int(key[2:])].name if key.startswith("ch") else key
        )

        if hit_tier_name == "pht":
            meta = detmeta.results.partition_ecal.cuspEmax_ctc_cal.eres_linear
        elif hit_tier_name == "hit":
            meta = detmeta.results.ecal.cuspEmax_ctc_cal.eres_linear

        out_dict[hpge] = meta

    return AttrsDict(out_dict)


def build_energy_res_func(function: str) -> Callable:
    """Energy resolution function builder."""
    if function == "FWHMLinear":
        return lambda energy, a, b: (a + b * energy) ** 0.5
    if function == "FWHMQuadratic":
        return lambda energy, a, b, c: (a + b * energy + c * energy * energy) ** 0.5

    raise NotImplementedError


def build_energy_res_func_dict(
    l200data: str | Path,
    metadata: LegendMetadata,
    runid: str,
    *,
    hit_tier_name: str = "hit",
    energy_res_pars: dict | AttrsDict | None = None,
) -> dict[str, Callable]:
    r"""Build energy resolution functions for each HPGe detector in a LEGEND-200 run.

    Returns
    -------
    Mapping of HPGe name to energy resolution function (FWHM), where energy is
    expected in units of keV.

    Parameters
    ----------
    l200data
        The path to the L200 data production cycle.
    metadata
        The metadata instance
    runid
        LEGEND-200 run identifier, must be of the form `{EXPERIMENT}-{PERIOD}-{RUN}-{TYPE}`.
    hit_tier_name
        name of the hit tier. This is typically "hit" or "pht".
    pars_db
        optional existing *non-lazy* instance of
        ``TextDB(".../path/to/prod/generated/par_{hit_tier_name}")``.
    """
    if energy_res_pars is None:
        energy_res_pars = lookup_energy_res_metadata(
            l200data,
            metadata,
            runid,
            hit_tier_name=hit_tier_name,
        )

    if not isinstance(energy_res_pars, AttrsDict):
        energy_res_pars = AttrsDict(energy_res_pars)

    _func_full = build_energy_res_func("FWHMLinear")

    energy_res_sigma_func = {}
    for hpge, meta in energy_res_pars.items():
        # use functools.partial correctly freeze the parameters into the function
        base = functools.partial(
            _func_full,
            a=meta.parameters.a,
            b=meta.parameters.b,
        )

        def _eres(E, base=base):
            return base(E)

        msg = f"measured FWHM for {hpge} at 2 MeV is ~{_eres(2000)} keV"
        log.debug(msg)

        energy_res_sigma_func[hpge] = _eres

    return energy_res_sigma_func
