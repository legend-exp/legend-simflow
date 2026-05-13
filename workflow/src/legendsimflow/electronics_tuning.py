"""Tune the electronics response parameters of the simulation against data superpulses.

Fits the Gaussian sigma and exponential tau of the system response kernel
by minimising the mean RMS between simulated and measured current
superpulses across drift-time slices.

"""

from __future__ import annotations

import logging

import awkward as ak
import numpy as np
from iminuit import Minuit
from scipy.interpolate import interp1d

from legendsimflow import psl
from legendsimflow.superpulses import Slice, Superpulse

log = logging.getLogger(__name__)

# Match make_hpge_realistic_pulse_shape_lib defaults
DT_DATA: float = 16.0
MW_PARS: dict[str, int] = {"length": 48, "num_mw": 3, "mw_type": 0}


def process_ideal_waveforms(  ### Merge this with psl.make_realistic_pulse_shape_lib to ensure consistency!!
    ideal_wfs: np.ndarray,
    rf_kernel: np.ndarray,
    dt: float,
    alignment_idx: int,
    nsamples_output: int,
) -> np.ndarray:
    """Apply the electronics chain to ideal charge waveforms.

    Mirrors the data DSP current estimation chain:
    convolve - downsample to data sampling (16 ns) - differentiate -
    upsample back to 1 ns - moving-window average - align to current peak.

    The downsample/upsample ensures the simulated current
    waveform has the same effective bandwidth as the data, where the
    charge waveform is acquired at 16 ns and the current is upsampled
    to 1 ns before smoothing.

    Parameters
    ----------
    ideal_wfs
        Ideal charge waveforms, shape ``(n_wfs, n_samples)``.
    rf_kernel
        System response kernel from
        :func:`psl.build_electronics_response_kernel`.
    dt
        Time step of the ideal waveforms in ns.
    alignment_idx
        Sample index where the current peak is placed after alignment.
    nsamples_output
        Length of the output waveforms.

    Returns
    -------
    aligned
        Aligned current waveforms, shape ``(n_wfs, nsamples_output)``.

    """
    from dspeed.processors import moving_window_multi, upsampler  # noqa: PLC0415

    convolved = ak.to_numpy(psl.apply_electronics_response(ideal_wfs, rf_kernel))

    # Downsample to data sampling rate (16 ns)
    downsample_factor = round(DT_DATA / dt)
    n_wfs, n_samples = convolved.shape
    # Trim to a length divisible by the downsample factor
    n_trim = (n_samples // downsample_factor) * downsample_factor
    convolved_trimmed = convolved[:, :n_trim]
    # Average bins
    charge_16ns = convolved_trimmed.reshape(n_wfs, -1, downsample_factor).mean(axis=-1)

    # Differentiate at 16 ns
    current_16ns = np.diff(charge_16ns, axis=-1, prepend=0)

    # Upsample x16 back to 1 ns
    n_samples_16 = current_16ns.shape[-1]
    n_samples_up = (n_samples_16 - 1) * downsample_factor
    current_1ns = np.zeros((n_wfs, n_samples_up), dtype=float)
    upsampler(
        current_16ns.astype(float, copy=False),
        downsample_factor,
        current_1ns,
    )

    # Moving window average
    mwa_out = np.zeros_like(current_1ns)
    moving_window_multi(
        current_1ns,
        MW_PARS["length"],
        MW_PARS["num_mw"],
        MW_PARS["mw_type"],
        mwa_out,
    )

    aligned, _ = psl.align_waveforms_to_peak(mwa_out, alignment_idx, nsamples_output)
    return aligned


def select_ideal_wfs_in_slice(
    ideal_wfs: np.ndarray, dt: float, sl: Slice
) -> np.ndarray:
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
    sim_avg: np.ndarray,
    sim_time: np.ndarray,
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
        return 1e6

    return float(np.sqrt(np.mean((sim_on_data - data_wf) ** 2)))


def build_cost_function(
    ideal_wfs_slice: dict[Slice, np.ndarray],
    data_superpulses: dict[Slice, Superpulse],
    dt: float,
    alignment_idx: int,
    nsamples_output: int,
    comparison_window: tuple[float, float] | None = None,
) -> callable:
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

        try:
            rf = psl.build_electronics_response_kernel(
                dt, mu_bandwidth=0.0, sigma_bandwidth=sigma, tau_rc=tau
            )
        except ValueError:
            return 1e6

        total = 0.0
        for sl, ideal_wfs in ideal_wfs_slice.items():
            processed = process_ideal_waveforms(
                ideal_wfs, rf, dt, alignment_idx, nsamples_output
            )
            sim_avg = np.mean(processed, axis=0)
            sim_time = (np.arange(len(sim_avg)) - alignment_idx) * dt
            total += compute_rms_in_slice(
                sim_avg, sim_time, data_superpulses[sl], comparison_window
            )
        return total / len(ideal_wfs_slice)

    return cost


def fit_electronics_parameters(
    ideal_pulse_shape_lib,
    data_superpulses: dict[Slice, Superpulse],
    angle: str = "000",
    *,
    sigma_start: float = 10.0,
    tau_start: float = 50.0,
    sigma_limits: tuple[float, float] = (0.1, 100),
    tau_limits: tuple[float, float] = (1, 500),
    comparison_window: tuple[float, float] | None = None,
    max_calls: int = 5000,
) -> dict:
    """Fit the electronics response parameters sigma and tau.

    Selects ideal waveforms from the pulse shape library in the same
    drift-time slices as the data superpulses, then minimises the mean
    RMS across slices using Minuit (MIGRAD + HESSE).

    Parameters
    ----------
    ideal_pulse_shape_lib
        Ideal waveform map (LGDO Struct) as read from LH5. Must contain
        ``waveform_{angle}_deg`` and ``dt``.
    data_superpulses
        Data superpulses keyed by slice, as returned by
        :func:`superpulses.read_superpulses_from_lh5`.
    angle
        Crystal axis angle tag, e.g. ``"000"``.
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
    run_minos
        If True, compute asymmetric profile-likelihood errors after HESSE.

    Returns
    -------
    result
        Dictionary with keys:

        - ``sigma``, ``tau`` : best-fit values in ns
        - ``best_rms`` : cost function value at the minimum
        - ``ideal_wfs_slice`` : ``dict[Slice, np.ndarray]`` of ideal
          charge waveforms selected in each drift-time slice
        - ``dt`` : time step of the ideal waveforms in ns
        - ``alignment_idx`` : sample index used for current-peak alignment
        - ``nsamples_output`` : length of the output current waveforms
          (derived from data superpulses; assumes all slices have the
          same current waveform length)
        - ``minuit`` : the :class:`iminuit.Minuit` object
        - ``cost_history`` : list of ``((sigma, tau), rms)`` per evaluation

    """
    from reboost import units  # noqa: PLC0415

    # Read library
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

    # Derive nsamples_output from data superpulses
    first_sp = next(iter(data_superpulses.values()))
    nsamples_output = len(first_sp.current_wf)

    # Slice ideal waveforms to match data superpulses
    ideal_wfs_slice: dict[Slice, np.ndarray] = {}
    for sl in data_superpulses:
        wfs = select_ideal_wfs_in_slice(ideal_wfs_flat, dt, sl)
        if len(wfs) > 0:
            ideal_wfs_slice[sl] = wfs

    if not ideal_wfs_slice:
        msg = "no valid slices found"
        raise RuntimeError(msg)

    log.info("prepared %d slices", len(ideal_wfs_slice))

    # Build cost
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

    # Minuit
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
    from matplotlib import pyplot as plt  # noqa: PLC0415

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
) -> tuple:
    """Overlay data and best-fit simulated current superpulses.

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

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : array of matplotlib.axes.Axes
    """
    from matplotlib import pyplot as plt  # noqa: PLC0415

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

    for idx, sl in enumerate(sorted_slices):
        ax = axes_flat[idx]
        data_sp = data_superpulses[sl]

        # Produce sim superpulse at best-fit parameters
        processed = process_ideal_waveforms(
            ideal_wfs_slice[sl], rf, dt, alignment_idx, nsamples_output
        )
        sim_avg = np.mean(processed, axis=0)
        sim_time = (np.arange(len(sim_avg)) - alignment_idx) * dt

        # Interpolate sim onto data grid for overlay and RMS
        f_interp = interp1d(
            sim_time, sim_avg, kind="linear", bounds_error=False, fill_value=np.nan
        )
        data_time = data_sp.current_time_axis
        data_wf = data_sp.current_wf

        ax.plot(data_time, data_wf, color="#2166ac", linewidth=1.8, label="Data")
        ax.plot(
            data_time,
            f_interp(data_time),
            color="#b2182b",
            linewidth=1.8,
            linestyle="--",
            label="Simulation",
        )

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
            f"dt = [{dt_lo:.0f}, {dt_hi:.0f}] ns\nRMS = {rms:.5f}", fontsize=10
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

    return fig, axes
