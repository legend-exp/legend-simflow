# Copyright (C) 2026 Giovanna Saleh <giovanna.saleh@studenti.unipd.it>
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
"""Tune the electronics response parameters of the simulation against data superpulses.

Fits the Gaussian sigma and exponential tau of the system response kernel
by minimising the mean RMS between simulated and measured current
superpulses across drift-time slices.
"""

from __future__ import annotations

import logging
from collections.abc import Callable

import numpy as np
from iminuit import Minuit
from lgdo import Struct
from matplotlib import pyplot as plt
from numpy.typing import NDArray
from reboost import units
from scipy.interpolate import interp1d

from legendsimflow import psl
from legendsimflow.superpulses import Slice, Superpulse

log = logging.getLogger(__name__)


def select_ideal_wfs_in_slice(ideal_wfs: NDArray, dt: float, sl: Slice) -> NDArray:
    """Select ideal waveforms whose drift time falls in a slice.

    Drift times are computed on the fly for the provided waveforms.

    Parameters
    ----------
    ideal_wfs
        Ideal charge waveforms, shape ``(n_wfs, n_samples)``.
    dt
        Time step in ns.
    sl
        Drift-time slice. Only ``sl.drift_time_range`` is used.

    Returns
    -------
    selected
        Waveforms in the slice, shape ``(n_selected, n_samples)``.

    """
    # Compute drift times
    raw_current = np.diff(ideal_wfs, axis=-1, prepend=0)
    nan_mask = np.all(np.isnan(raw_current), axis=-1)
    peak_idx = np.argmax(np.nan_to_num(raw_current, nan=0.0), axis=-1)
    drift_times = np.where(nan_mask, np.nan, peak_idx * dt)

    # Select waveforms in slice
    lo, hi = sl.drift_time_range
    mask = (drift_times >= lo) & (drift_times < hi)
    selected = ideal_wfs[mask]

    if len(selected) == 0:
        log.warning("no ideal waveforms in %s", sl)

    return selected


def compute_rms_in_slice(
    sim_avg: NDArray,
    sim_time: NDArray,
    data_sp: Superpulse,
    comparison_window: tuple[float, float] | None = None,
) -> float:
    """RMS residual between a simulated and a data current superpulse.

    The simulation is linearly interpolated onto the data time grid.
    Samples outside the sim time range are excluded from the comparison.

    Parameters
    ----------
    sim_avg
        Simulated average current waveform, shape ``(n_sim,)``.
    sim_time
        Time axis of the simulation in ns, shape ``(n_sim,)``.
    data_sp
        Data superpulse. Only ``current_wf`` and ``current_time_axis``
        are read.
    comparison_window
        ``(t_min, t_max)`` in ns relative to the current peak.
        If ``None``, the full waveform overlap is used.

    Returns
    -------
    rms
        Root mean square of the residuals.

    """
    data_time = data_sp.current_time_axis
    data_wf = data_sp.current_wf

    f = interp1d(
        sim_time, sim_avg, kind="linear", bounds_error=False, fill_value=np.nan
    )
    sim_on_data = f(data_time)

    if comparison_window is not None:
        lo, hi = comparison_window
        mask = (data_time >= lo) & (data_time <= hi)
        sim_on_data = sim_on_data[mask]
        data_wf = data_wf[mask]

    valid = np.isfinite(sim_on_data)
    sim_on_data = sim_on_data[valid]
    data_wf = data_wf[valid]

    if len(data_wf) == 0:
        msg = "no valid samples in comparison region"
        raise ValueError(msg)

    return float(np.sqrt(np.mean((sim_on_data - data_wf) ** 2)))


def build_cost_function(
    ideal_wfs_slice: dict[Slice, NDArray],
    data_superpulses: dict[Slice, Superpulse],
    dt: float,
    alignment_idx: int,
    nsamples_output: int,
    comparison_window: tuple[float, float] | None = None,
) -> Callable:
    """Build the scalar cost function for the Minuit minimiser.

    The returned function has signature ``cost(sigma, tau) -> float``
    and computes the mean RMS across all drift-time slices.

    Parameters
    ----------
    ideal_wfs_slice
        Preselected ideal waveforms per slice,
        ``{Slice: ideal_wfs_array}``.
    data_superpulses
        Data superpulses keyed by the same slices.
    dt
        Time step of the ideal waveforms in ns.
    alignment_idx
        Sample index where the current peak is placed after alignment.
    nsamples_output
        Length of the output waveforms.
    comparison_window
        ``(t_min, t_max)`` in ns. Passed through to :func:`compute_rms_in_slice`.

    Returns
    -------
    cost
        ``cost(sigma, tau) -> float``.

    """

    def cost(sigma, tau):
        if sigma <= 0 or tau <= 0:
            return 1e6

        rf = psl.build_electronics_response_kernel(
            dt, mu_bandwidth=0.0, sigma_bandwidth=sigma, tau_rc=tau
        )

        total = 0.0
        for sl, ideal_wfs in ideal_wfs_slice.items():
            processed, _ = psl.process_ideal_waveforms(
                ideal_wfs,
                rf,
                dt,
                alignment_idx,
                nsamples_output,
                mw_pars=psl.MW_PARS,
                dt_data=psl.DT_DATA,
            )
            sim_avg = np.mean(processed, axis=0)
            sim_time = (np.arange(len(sim_avg)) - alignment_idx) * dt
            total += compute_rms_in_slice(
                sim_avg, sim_time, data_superpulses[sl], comparison_window
            )
        return total / len(ideal_wfs_slice)

    return cost


def get_ideal_wfs_all_slices(
    ideal_pulse_shape_lib: Struct,
    data_superpulses: dict[Slice, Superpulse],
    angle: str = "000",
) -> dict:
    """Select ideal waveforms per drift-time slice.

    Reads the ideal pulse shape library, flattens the (r, z) grid,
    and selects waveforms whose drift time falls in each data
    superpulse slice.

    Parameters
    ----------
    ideal_pulse_shape_lib
        Ideal waveform map (LGDO Struct) as read from lh5. Must
        contain ``waveform_{angle}_deg`` (an Array) and ``dt`` (time step), a Scalar.
    data_superpulses
        Data superpulses keyed by slice.
    angle
        Crystal axis angle tag, e.g. ``"000"``.

    Returns
    -------
    dict
        Keys:

        - ``ideal_wfs_slice`` : ``dict[Slice, NDArray]``
        - ``dt`` : time step in ns
        - ``alignment_idx`` : sample index for current-peak alignment
        - ``nsamples_output`` : output waveform length (from data)

    """
    wf_key = f"waveform_{angle}_deg"
    ideal_wfs_full = ideal_pulse_shape_lib[wf_key].view_as("np")
    dt = ideal_pulse_shape_lib["dt"].value * units.units_convfact(
        ideal_pulse_shape_lib["dt"], "ns"
    )

    # Flatten (n_r, n_z, n_samples) -> (n_pixels, n_samples), drop NaN pixels
    orig_shape = ideal_wfs_full.shape
    ideal_wfs_flat = ideal_wfs_full.reshape(-1, orig_shape[-1])
    valid_mask = ~np.any(np.isnan(ideal_wfs_flat), axis=-1)
    ideal_wfs_flat = ideal_wfs_flat[valid_mask]

    alignment_idx = orig_shape[-1] // 2

    first_sp = next(iter(data_superpulses.values()))
    nsamples_output = len(first_sp.current_wf)

    ideal_wfs_slice: dict[Slice, NDArray] = {}
    for sl in data_superpulses:
        wfs = select_ideal_wfs_in_slice(ideal_wfs_flat, dt, sl)
        if len(wfs) > 0:
            ideal_wfs_slice[sl] = wfs

    if not ideal_wfs_slice:
        msg = "no valid slices found"
        raise RuntimeError(msg)

    log.info("prepared %d slices", len(ideal_wfs_slice))

    return {
        "ideal_wfs_slice": ideal_wfs_slice,
        "dt": dt,
        "alignment_idx": alignment_idx,
        "nsamples_output": nsamples_output,
    }


def fit_electronics_parameters(
    ideal_wfs_slice: dict[Slice, NDArray],
    data_superpulses: dict[Slice, Superpulse],
    dt: float,
    alignment_idx: int,
    nsamples_output: int,
    *,
    sigma_start: float,
    tau_start: float,
    sigma_limits: tuple[float, float],
    tau_limits: tuple[float, float],
    comparison_window: tuple[float, float] | None = None,
    max_calls: int = 5000,
) -> dict:
    """Fit the electronics response parameters sigma and tau.

    Minimises the mean RMS between simulated and measured current
    superpulses across drift-time slices using Minuit (MIGRAD).

    Parameters
    ----------
    ideal_wfs_slice
        Ideal charge waveforms per slice, as returned by
        :func:`get_ideal_wfs_all_slices`.
    data_superpulses
        Data superpulses keyed by the same slices.
    dt
        Time step of the ideal waveforms in ns.
    alignment_idx
        Sample index for current-peak alignment.
    nsamples_output
        Length of the output waveforms.
    sigma_start
        Initial value for the Gaussian sigma in ns.
    tau_start
        Initial value for the exponential tau in ns.
    sigma_limits
        Hard bounds ``(lo, hi)`` for sigma in ns.
    tau_limits
        Hard bounds ``(lo, hi)`` for tau in ns.
    comparison_window
        ``(t_min, t_max)`` in ns relative to the current peak.
        If ``None``, the full waveform overlap is used.
    max_calls
        Maximum number of Minuit function evaluations.

    Returns
    -------
    dict
        Keys: ``sigma``, ``tau``, ``best_rms``,
        ``ideal_wfs_slice``, ``dt``, ``alignment_idx``,
        ``nsamples_output``, ``minuit``, ``cost_history``.

    """
    cost_fn = build_cost_function(
        ideal_wfs_slice,
        data_superpulses,
        dt,
        alignment_idx,
        nsamples_output,
        comparison_window,
    )

    history: list[tuple[tuple[float, float], float]] = []

    def tracked_cost(sigma, tau):
        val = cost_fn(sigma, tau)
        history.append(((sigma, tau), val))
        return val

    m = Minuit(tracked_cost, sigma=sigma_start, tau=tau_start)
    m.errors = (5.0, 10.0)
    m.limits["sigma"] = sigma_limits
    m.limits["tau"] = tau_limits

    m.migrad(ncall=max_calls)
    if not m.valid:
        log.warning("MIGRAD did not converge")

    return {
        "sigma": m.values["sigma"],
        "tau": m.values["tau"],
        "best_rms": m.fval,
        "ideal_wfs_slice": ideal_wfs_slice,
        "dt": dt,
        "alignment_idx": alignment_idx,
        "nsamples_output": nsamples_output,
        "minuit": m,
        "cost_history": history,
    }


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


def plot_convergence(result: dict) -> tuple:
    """Convergence diagnostics for the Minuit optimisation.

    Three panels: RMS vs function call, sigma trajectory, tau trajectory.
    The best-fit value is marked with a dashed red line in each parameter
    panel.

    Parameters
    ----------
    result
        Dictionary returned by :func:`fit_electronics_parameters`.
        Uses ``cost_history``, ``sigma``, ``tau``, and ``best_rms``.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : tuple of matplotlib.axes.Axes
    """
    history = result["cost_history"]
    calls = np.arange(len(history))
    params = np.array([h[0] for h in history])  # (N, 2) - sigma, tau
    rms = np.array([h[1] for h in history])

    # Clip outliers for display (penalty returns 1e6)
    valid = rms < 1e5
    rms_display = np.where(rms < 100000.0, rms, np.nan) if valid.any() else rms

    fig, (ax_rms, ax_sigma, ax_tau) = plt.subplots(
        1, 3, figsize=(15, 4.5), constrained_layout=True
    )

    # RMS convergence
    ax_rms.plot(calls, rms_display, linewidth=0.8, color="#1878b2", alpha=0.8)
    ax_rms.axhline(
        result["best_rms"],
        color="#b2182b",
        linewidth=1.2,
        linestyle="--",
        label=f"best = {result['best_rms']:.5f}",
    )
    ax_rms.set_xlabel("Function call")
    ax_rms.set_ylabel("Mean RMS")
    ax_rms.set_title("Cost convergence")
    ax_rms.set_yscale("log")
    ax_rms.legend(fontsize=9)
    ax_rms.grid(visible=True, linestyle="--", alpha=0.3)

    # Sigma
    ax_sigma.plot(calls, params[:, 0], linewidth=0.8, color="#e08214", alpha=0.8)
    ax_sigma.axhline(
        result["sigma"],
        color="#b2182b",
        linewidth=1.2,
        linestyle="--",
        label=rf"best = {result['sigma']:.2f} ns",
    )
    ax_sigma.set_xlabel("Function call")
    ax_sigma.set_ylabel(r"$\sigma$ [ns]")
    ax_sigma.set_title(r"$\sigma$ trajectory")
    ax_sigma.legend(fontsize=9)
    ax_sigma.grid(visible=True, linestyle="--", alpha=0.3)

    # Tau
    ax_tau.plot(calls, params[:, 1], linewidth=0.8, color="#542788", alpha=0.8)
    ax_tau.axhline(
        result["tau"],
        color="#b2182b",
        linewidth=1.2,
        linestyle="--",
        label=rf"best = {result['tau']:.2f} ns",
    )
    ax_tau.set_xlabel("Function call")
    ax_tau.set_ylabel(r"$\tau$ [ns]")
    ax_tau.set_title(r"$\tau$ trajectory")
    ax_tau.legend(fontsize=9)
    ax_tau.grid(visible=True, linestyle="--", alpha=0.3)

    return fig, (ax_rms, ax_sigma, ax_tau)


def plot_best_fit(
    result: dict,
    data_superpulses: dict[Slice, Superpulse],
    comparison_window: tuple[float, float] | None = None,
    plot_window: tuple[float, float] | None = None,
    plot_charge: bool = False,
    detector: str| None = None
) -> tuple:
    """Overlay data and best-fit simulated superpulses (current or charge).

    One panel per drift-time slice, sorted by drift time. Each panel
    shows the data superpulse and the simulation at the best-fit
    (sigma, tau), with the per-slice RMS annotated.

    Parameters
    ----------
    result
        Dictionary returned by :func:`fit_electronics_parameters`.
        Uses ``sigma``, ``tau``, ``best_rms``, ``ideal_wfs_slice``,
        ``dt``, ``alignment_idx``, and ``nsamples_output``.
    data_superpulses
        Data superpulses keyed by slice.
    comparison_window
        ``(t_min, t_max)`` in ns relative to the current peak. If given,
        the window is shaded on each panel.
    plot_window
        ``(t_min, t_max)`` in ns for the x-axis limits.
        Defaults to ``comparison_window`` if set, otherwise auto-scaled.
    plot_charge
        Plot the charge (instead of the current), waveforms.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : array of matplotlib.axes.Axes
    data_amax: the Amax value of the highest drift time slice
    mc_amax: the Amax value of the highest drift time slice
    """
    sigma = result["sigma"]
    tau = result["tau"]
    ideal_wfs_slice = result["ideal_wfs_slice"]
    dt = result["dt"]
    alignment_idx = result["alignment_idx"]
    nsamples_output = result["nsamples_output"]

    rf = psl.build_electronics_response_kernel(
        dt, mu_bandwidth=0.0, sigma_bandwidth=sigma, tau_rc=tau
    )

    # Sort slices by drift-time centre
    sorted_slices = sorted(ideal_wfs_slice.keys(), key=lambda s: s.drift_time_center)
    n_slices = len(sorted_slices)

    ncols = min(4, n_slices)
    nrows = int(np.ceil(n_slices / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.5 * ncols, 4 * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    axes_flat = axes.flatten()
    data_amax = None
    mc_amax = None

    for idx, sl in enumerate(sorted_slices):
        ax = axes_flat[idx]
        data_sp = data_superpulses[sl]

        # Produce sim superpulse at best-fit parameters
        processed, _ = psl.process_ideal_waveforms(
            ideal_wfs_slice[sl],
            rf,
            dt,
            alignment_idx,
            nsamples_output,
            mw_pars=psl.MW_PARS,
            dt_data=psl.DT_DATA,
            return_mode="charge" if plot_charge else "current",
        )

        sim_avg = np.mean(processed, axis=0)
        sim_time = (np.arange(len(sim_avg)) - alignment_idx) * dt

        # Interpolate sim onto data grid for overlay and RMS
        f_interp = interp1d(
            sim_time, sim_avg, kind="linear", bounds_error=False, fill_value=np.nan
        )

        data_time = (
            data_sp.current_time_axis if not plot_charge else data_sp.charge_time_axis
        )
        data_wf = data_sp.current_wf if not plot_charge else data_sp.charge_wf

        ax.plot(data_time, data_wf, color="#2166ac", linewidth=1.8, label="Data")
        ax.plot(
            data_time,
            f_interp(data_time),
            color="#b2182b",
            linewidth=1.8,
            linestyle="--",
            label="Simulation",
        )
        if (idx == len(sorted_slices) - 1) and not plot_charge:
            data_amax = float(np.max(data_wf))
            mc_amax = float(np.max(f_interp(data_time)))

        # Shade comparison window
        if comparison_window is not None:
            ax.axvspan(
                comparison_window[0],
                comparison_window[1],
                alpha=0.08,
                color="grey",
            )

        # Per-slice RMS
        rms = compute_rms_in_slice(sim_avg, sim_time, data_sp, comparison_window)
        dt_lo, dt_hi = sl.drift_time_range
        ax.set_title(
            f"{detector} dt = [{dt_lo:.0f}, {dt_hi:.0f}] ns\nRMS = {rms:.5f}", fontsize=10
        )

        # Axis limits
        if plot_window is not None:
            ax.set_xlim(plot_window)
        elif comparison_window is not None:
            ax.set_xlim(comparison_window)

        ax.set_xlabel("Time [ns]")
        ax.grid(visible=True, linestyle="--", alpha=0.3)

        if idx % ncols == 0:
            ax.set_ylabel("Amplitude [a.u.]")
        if idx == 0:
            ax.legend(fontsize=9)

    # Hide unused panels
    for ax in axes_flat[n_slices:]:
        ax.set_visible(False)

    fig.suptitle(
        rf"Best fit: $\sigma$ = {sigma:.2f} ns, $\tau$ = {tau:.2f} ns"
        f"  |  mean RMS = {result['best_rms']:.5f}",
        fontsize=12,
    )

    return fig, axes, data_amax, mc_amax
