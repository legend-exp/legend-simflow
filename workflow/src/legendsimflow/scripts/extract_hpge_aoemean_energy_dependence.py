# ruff: noqa: I002

# Copyright (C) 2026 Toby Dixon <toby.dixon.23@ucl.ac.uk>
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

"""Extract the HPGe A/E mean energy dependence from electron-gun simulations.

For every modelable HPGe detector of a run, the A/E of mono-energetic electrons
uniformly generated in the bulk of the detectors (one simulation per energy, see
:mod:`.simulate_electron_gun`) is computed with the same PSD routines used by
the ``hit`` tier, without applying any energy-dependence correction. Electrons
are single-site events, so the raw A/E distribution at each energy is peaked;
its position is estimated with the unbinned half-sample mode, which ignores the
low-side tail left by bremsstrahlung losses. The modes are fitted as a function
of the energy with the same linear model and fitting procedure used for
LEGEND-200 data
(:class:`pygama.pargen.AoE_cal.Pol1`), which is the A/E mean model consumed by
the ``hit`` tier (see :func:`legendsimflow.hpge_pars.build_aoe_mean_func_dict`).
"""

import argparse
import logging
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
import lh5
import numpy as np
import pyg4ometry
import pygeomhpges
import pygeomtools
from dbetto import AttrsDict
from dbetto.utils import load_dict
from iminuit import Minuit, cost
from lh5 import LH5Iterator
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from numpy.typing import NDArray
from pygama.pargen.AoE_cal import Pol1
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import metadata as mutils
from legendsimflow import nersc, patterns, utils
from legendsimflow import reboost as reboost_utils
from legendsimflow.metadata import get_tier_settings
from legendsimflow.plot import decorate
from legendsimflow.scripts import log_script_invocation

# style of the PSD methods in the validation plots
_SIM_TYPE_STYLE = {
    "single_template": {"color": "#0077BB", "label": "single-temp"},
    "psl": {"color": "#CC3311", "label": "pulse lib"},
}

# scale factor from the median absolute deviation to the standard deviation of
# a Gaussian sample
_MAD_TO_SIGMA = 1.4826

# number of bootstrap replicas used to estimate the uncertainty on the mode
_N_BOOTSTRAP = 200

# the energy at which the intrinsic resolution is reported (QA)
_RESOLUTION_REF_ENERGY_KEV = 2000


def group_stp_files_by_energy(files: Iterable[str | Path]) -> dict[int, list[Path]]:
    """Group the electron-gun ``stp`` files by electron energy.

    The energy is read from the file name (see
    :func:`legendsimflow.patterns.electron_gun_energy_from_path`).

    Returns
    -------
    dict[int, list[Path]]
        Mapping of electron energy (keV) to the sorted list of files, sorted by
        energy.
    """
    out: dict[int, list[Path]] = {}
    for f in files:
        path = Path(f)
        out.setdefault(patterns.electron_gun_energy_from_path(path), []).append(path)
    return {e: sorted(out[e]) for e in sorted(out)}


def half_sample_mode(sample: NDArray) -> float:
    """Unbinned mode of a sample, with the half-sample mode estimator.

    Recursively selects the half of the (sorted) sample with the smallest
    range, until three or fewer values are left, and returns their centre. The
    estimator needs no binning and no bandwidth, and is insensitive to the
    tails of the distribution, which is what makes it suitable for the raw A/E
    of the electron-gun simulations: a peak with a low-side tail from
    bremsstrahlung losses.

    Parameters
    ----------
    sample
        The values whose mode is estimated. Must be non-empty and finite.

    Returns
    -------
    float
        The estimated mode.
    """
    x = np.sort(np.asarray(sample, dtype=float))
    n = x.size
    if n == 0:
        msg = "cannot estimate the mode of an empty sample"
        raise ValueError(msg)

    while n > 3:
        # smallest interval containing half of the remaining values
        m = (n + 1) // 2
        widths = x[m - 1 :] - x[: n - m + 1]
        i = int(np.argmin(widths))
        x = x[i : i + m]
        n = m

    if n == 1:
        return float(x[0])
    if n == 2:
        return float(x.mean())
    # three values left: the mode is the centre of the closest pair
    low, high = x[1] - x[0], x[2] - x[1]
    if low < high:
        return float((x[0] + x[1]) / 2)
    if high < low:
        return float((x[1] + x[2]) / 2)
    return float(x[1])


def summarize_aoe(aoe: NDArray, rng: np.random.Generator | None = None) -> dict:
    """Summary statistics of a raw A/E sample at fixed electron energy.

    Non-finite values are discarded. The position of the distribution is
    estimated with :func:`half_sample_mode`, and its uncertainty with a
    bootstrap over ``200`` replicas (the estimator has no simple analytic
    variance). The median, the standard deviation and the median absolute
    deviation are also reported, for quality assurance.

    Parameters
    ----------
    aoe
        The raw A/E values.
    rng
        Random generator used by the bootstrap. Seeded deterministically if
        not given, so that the extraction is reproducible.

    Returns
    -------
    dict
        Keys ``n_events``, ``mode``, ``mode_err``, ``median``, ``std``, ``mad``.
    """
    aoe = np.asarray(aoe, dtype=float)
    aoe = aoe[np.isfinite(aoe)]
    n = int(aoe.size)
    if n == 0:
        return dict.fromkeys(("mode", "mode_err", "median", "std", "mad"), np.nan) | {
            "n_events": 0
        }

    if rng is None:
        rng = np.random.default_rng(1234)

    mode = half_sample_mode(aoe)

    # the sampling distribution of the half-sample mode has no simple closed
    # form: estimate its width by resampling
    replicas = [
        half_sample_mode(rng.choice(aoe, size=n, replace=True))
        for _ in range(_N_BOOTSTRAP)
    ]
    mode_err = float(np.std(replicas, ddof=1))

    median = float(np.median(aoe))
    return {
        "n_events": n,
        "mode": mode,
        "mode_err": mode_err,
        "median": median,
        "std": float(np.std(aoe)),
        "mad": float(np.median(np.abs(aoe - median))),
    }


def fit_linear_energy_dependence(
    energies: Sequence[float], values: Sequence[float], errors: Sequence[float]
) -> tuple[dict[str, float], dict[str, float]]:
    """Fit the A/E energy dependence with the linear model used on data.

    Mirrors the mean fit of
    :meth:`pygama.pargen.AoE_cal.CalAoE.energy_correction`: the
    :class:`pygama.pargen.AoE_cal.Pol1` model (``x*a+b``) and its initial
    guess, a least-squares cost with the soft-L1 loss, and the simplex, migrad
    and hesse Minuit steps. Points with non-finite value or error are
    discarded; a floor is applied to the errors to keep the weights finite when
    the sample is degenerate.

    Returns
    -------
    pars, errs
        Best-fit parameters and their uncertainties, keyed by ``a`` and ``b``.

    Raises
    ------
    ValueError
        If fewer than two valid points are available.
    RuntimeError
        If the fit does not converge.
    """
    x = np.asarray(energies, dtype=float)
    y = np.asarray(values, dtype=float)
    ey = np.asarray(errors, dtype=float)

    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(ey)
    x, y, ey = x[ok], y[ok], ey[ok]
    if x.size < 2:
        msg = f"at least two valid points are needed for the linear fit, got {x.size}"
        raise ValueError(msg)

    ey = np.maximum(ey, 1e-9)

    c = cost.LeastSquares(x, y, ey, Pol1.func)
    c.loss = "soft_l1"
    m = Minuit(c, *Pol1.guess(x, y, ey))
    m.simplex()
    m.migrad()
    m.hesse()

    if not m.valid:
        msg = "the linear fit of the A/E modes vs. energy did not converge"
        raise RuntimeError(msg)

    pars = {"a": float(m.values["a"]), "b": float(m.values["b"])}
    errs = {"a": float(m.errors["a"]), "b": float(m.errors["b"])}
    return pars, errs


def extract_energy_dependence(
    samples: Mapping[int, NDArray], log_prefix: str, log: logging.Logger
) -> tuple[
    dict | None, dict[str, list], dict[int, dict[str, float | int]], dict | None
]:
    """Summarise the per-energy A/E samples of one detector and PSD method and fit them.

    Parameters
    ----------
    samples
        Mapping of electron energy (keV) to the raw A/E sample.
    log_prefix
        Prefix of the log messages (detector, run, method).
    log
        Logger.

    Returns
    -------
    model, points, stats, resolution
        The fitted ``{expression, pars, errs}`` model (``None`` if the fit is
        not possible), the per-energy statistics as arrays (keyed by
        statistic) and as a mapping keyed by energy, and the intrinsic A/E
        resolution (in percent) at the energy closest to 2 MeV.
    """
    energies = sorted(samples)
    stats = {e: summarize_aoe(samples[e]) for e in energies}
    points = {
        "energy_in_keV": [int(e) for e in energies],
        **{
            key: [stats[e][key] for e in energies]
            for key in ("n_events", "mode", "mode_err", "median", "std", "mad")
        },
    }
    for e in energies:
        log.info(
            "%s at %d keV: n = %d, mode = %.5f +- %.5f, std = %.2e",
            log_prefix,
            e,
            stats[e]["n_events"],
            stats[e]["mode"],
            stats[e]["mode_err"],
            stats[e]["std"],
        )

    model = None
    try:
        pars, errs = fit_linear_energy_dependence(
            points["energy_in_keV"], points["mode"], points["mode_err"]
        )
        model = {"expression": Pol1.string_func("x"), "pars": pars, "errs": errs}
        log.info(
            "%s: a = %.3e +- %.1e 1/keV, b = %.5f +- %.1e",
            log_prefix,
            pars["a"],
            errs["a"],
            pars["b"],
            errs["b"],
        )
    except (ValueError, RuntimeError) as e:
        log.warning("%s: could not fit the A/E energy dependence: %s", log_prefix, e)

    resolution = None
    valid = [e for e in energies if np.isfinite(stats[e]["std"])]
    if valid:
        e_ref = min(valid, key=lambda e: abs(e - _RESOLUTION_REF_ENERGY_KEV))
        resolution = {
            "val": 100 * float(stats[e_ref]["std"]),
            "energy_in_keV": int(e_ref),
        }

    return model, points, stats, resolution


def _plot_modes(
    pdf: PdfPages,
    points: Mapping[str, Mapping[str, list]],
    models: Mapping[str, Mapping],
    title: str,
) -> None:
    """Page with the A/E modes vs. energy and the linear fits."""
    fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
    e_line = np.linspace(
        min(mutils.ELECTRON_GUN_ENERGIES_IN_KEV) - 100,
        max(mutils.ELECTRON_GUN_ENERGIES_IN_KEV) + 100,
        200,
    )
    for sim_type, pts in points.items():
        style = _SIM_TYPE_STYLE[sim_type]
        ax.errorbar(
            pts["energy_in_keV"],
            pts["mode"],
            yerr=pts["mode_err"],
            marker="o",
            linestyle="none",
            capsize=2,
            color=style["color"],
            label=f"electron gun ({style['label']})",
        )
        if sim_type in models:
            pars = models[sim_type]["pars"]
            ax.plot(
                e_line,
                Pol1.func(e_line, **pars),
                color=style["color"],
                label=f"best fit ({style['label']}), a = {pars['a'] * 1e5:.3f} %/MeV",
            )
    ax.set_xlabel("electron energy [keV]")
    ax.set_ylabel("raw A/E mode")
    ax.legend()
    ax.set_title(title)
    decorate(fig)
    pdf.savefig()
    plt.close(fig)


def _plot_distributions(
    pdf: PdfPages,
    samples: Mapping[int, NDArray],
    stats: Mapping[int, Mapping[str, float]],
    title: str,
) -> None:
    """Page with the raw A/E distribution at each electron energy."""
    energies = sorted(samples)
    ncols = 3
    nrows = max(1, int(np.ceil(len(energies) / ncols)))
    fig, axs = plt.subplots(
        nrows,
        ncols,
        figsize=(4 * ncols, 3 * nrows),
        squeeze=False,
        layout="constrained",
    )
    for ax, energy in zip(axs.flat, energies, strict=False):
        aoe = np.asarray(samples[energy], dtype=float)
        aoe = aoe[np.isfinite(aoe)]
        st = stats[energy]
        if aoe.size > 0 and np.isfinite(st["mode"]):
            # robust width: the tails inflate the standard deviation
            width = 8 * max(_MAD_TO_SIGMA * st["mad"], 1e-6)
            ax.hist(
                aoe,
                bins=100,
                range=(st["mode"] - width, st["mode"] + width),
                histtype="stepfilled",
                alpha=0.6,
            )
            ax.axvline(st["mode"], color="black", linestyle="--", label="mode")
            ax.axvline(st["median"], color="grey", linestyle=":", label="median")
            ax.legend(
                title=f"n = {st['n_events']}\nmode = {st['mode']:.5f}\n"
                f"MAD = {st['mad']:.2e}",
                fontsize="small",
                title_fontsize="small",
            )
        ax.set_title(f"{energy} keV electrons")
        ax.set_xlabel("raw A/E")
        ax.set_ylabel("counts")
    for ax in list(axs.flat)[len(energies) :]:
        ax.set_axis_off()
    fig.suptitle(title)
    decorate(fig)
    pdf.savefig()
    plt.close(fig)


@snakemake_compatible(
    mapping={
        "runid": "wildcards.runid",
        "electron_stp_files": "input.electron_stp_files",
        "geom_file": "input.geom",
        "is_modelable_file": "input.is_modelable",
        "pars_file": "output.pars_file",
        "stats_file": "output.stats_file",
        "plot_file": "output.plot_file",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Extract the A/E mean energy dependence of the HPGe detectors of a "
            "run from electron-gun simulations."
        )
    )
    parser.add_argument("--runid", required=True, help="LEGEND run identifier")
    parser.add_argument(
        "--electron-stp-files",
        nargs="+",
        required=True,
        help=(
            "electron-gun stp tier files (all energies and jobs); the electron "
            "energy is read from the file name"
        ),
    )
    parser.add_argument("--geom-file", required=True, help="GDML geometry file")
    parser.add_argument(
        "--is-modelable-file",
        required=True,
        help="detinfo YAML file with the HPGe modeling status (runid -> detector -> flag)",
    )
    parser.add_argument(
        "--pars-file",
        required=True,
        help="output YAML file with the A/E mean energy-dependence models (per detector)",
    )
    parser.add_argument(
        "--stats-file",
        required=True,
        help="output YAML file with the per-detector electron-gun A/E statistics",
    )
    parser.add_argument("--plot-file", required=True, help="output PDF plot file")
    parser.add_argument("--log-file", default=None, help="log file")
    parser.add_argument(
        "--simflow-config",
        "--config",
        dest="simflow_config",
        required=True,
        help="simflow config YAML path",
    )

    args = parser.parse_args()

    config = utils.init_simflow_context(args.simflow_config, workflow=None).config
    metadata = config.metadata

    log = ldfs.utils.build_log(metadata.simprod.config.logging, args.log_file)
    log_script_invocation(log, "extract-hpge-aoemeanmod", parser, args)

    runid = args.runid
    gdml_file = nersc.dvs_ro(config, args.geom_file)
    stp_files = [nersc.dvs_ro(config, Path(f)) for f in args.electron_stp_files]

    tier_hit_settings = get_tier_settings(config, "hit")
    dead_layer_fraction = tier_hit_settings.dead_layer_fraction
    simulate_psd = tier_hit_settings.get("simulate_psd", True)
    simulate_psd_with_psl = tier_hit_settings.get("simulate_psd_with_psl", False)
    buffer_len = tier_hit_settings.buffer_len

    # the modelable detectors of this run
    is_modelable = load_dict(nersc.dvs_ro(config, args.is_modelable_file))
    dets = sorted(d for d, flag in is_modelable.get(runid, {}).items() if flag)
    log.info("modelable HPGe detectors in %s: %s", runid, ", ".join(dets))

    log.debug("loading the geometry")
    geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
    sens_tables = pygeomtools.detectors.get_all_senstables(geom)

    files_by_energy = group_stp_files_by_energy(stp_files)
    log.info(
        "electron energies found: %s keV", ", ".join(str(e) for e in files_by_energy)
    )
    # the detector origins are the same in all files
    det_loc_all = reboost_utils.read_detector_origins(
        next(iter(files_by_energy.values()))[0]
    )

    # no correction, no classifier: only the raw A/E is of interest here
    identity_psdcuts = AttrsDict({"aoe": {"low_side": 0.0, "high_side": 0.0}})

    # merged current-pulse models of the run, shared by all the detectors
    currmod_pars_all = (
        AttrsDict(
            load_dict(
                nersc.dvs_ro(
                    config, patterns.output_currmod_merged_filename(config, runid=runid)
                )
            )
        )
        if simulate_psd
        else AttrsDict({})
    )

    models: dict[str, dict] = {}
    stats_out: dict[str, dict] = {}

    plot_file = Path(args.plot_file)
    plot_file.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(str(plot_file)) as pdf:
        for det_idx, det in enumerate(dets):
            log.info("processing detector %s [%d/%d]", det, det_idx + 1, len(dets))

            if det not in sens_tables or sens_tables[det].detector_type != "germanium":
                log.warning(
                    "%s is not a germanium detector of the geometry, skipping", det
                )
                continue

            pyobj = pygeomhpges.make_hpge(
                sens_tables[det].metadata,
                registry=None,
                allow_cylindrical_asymmetry=False,
            )
            fccd = mutils.get_sanitized_fccd(metadata, det)
            det_loc = det_loc_all[det]

            # PSD inputs of this detector in this run, the same used by the hit
            # tier. free the previous detector's inputs before loading the next
            psd_inputs = None
            psd_inputs = reboost_utils.load_hpge_psd_inputs(
                config,
                det,
                runid,
                simulate_psd=simulate_psd,
                simulate_psd_with_psl=simulate_psd_with_psl,
                currmod_pars_all=currmod_pars_all,
            )
            methods = []
            if psd_inputs.can_model_psd:
                methods.append("single_template")
            if psd_inputs.can_model_psd_with_psl:
                methods.append("psl")
            if not methods:
                log.warning(
                    "the PSD response of %s in %s cannot be simulated (drift-time "
                    "map, current model or pulse-shape library missing), skipping",
                    det,
                    runid,
                )
                continue

            # collect the raw A/E samples per PSD method and electron energy
            samples: dict[str, dict[int, list[NDArray]]] = {
                st: {e: [] for e in files_by_energy} for st in methods
            }
            stp_table_name = f"stp/{det}"
            for energy, files in files_by_energy.items():
                # a job might have recorded no hit at all in this detector
                files_with_det = [f for f in files if lh5.ls(f, stp_table_name)]
                if len(files_with_det) == 0:
                    log.warning("no %s hits found in the %d keV files", det, energy)
                    continue

                iterator = LH5Iterator(
                    [str(f) for f in files_with_det],
                    stp_table_name,
                    buffer_len=buffer_len,
                )
                for lgdo_chunk in iterator:
                    chunk = lgdo_chunk.view_as("ak", with_units=True)

                    edep_active, energy_true = reboost_utils.hpge_active_energy(
                        chunk, pyobj, det_loc, fccd, dead_layer_fraction
                    )
                    has_energy = np.asarray(energy_true > 0)

                    psd, psd_detailed = reboost_utils.compute_hpge_psd_observables(
                        chunk,
                        edep_active,
                        energy_true,
                        det_loc,
                        psd_inputs,
                        aoe_res=1.0,
                        aoe_mean=1.0,
                        aoe_mean_psl=1.0,
                        psdcuts=identity_psdcuts,
                    )
                    if "single_template" in samples:
                        samples["single_template"][energy].append(
                            np.asarray(psd.aoe, dtype=float)[has_energy]
                        )
                    if "psl" in samples:
                        samples["psl"][energy].append(
                            np.asarray(psd_detailed.aoe, dtype=float)[has_energy]
                        )

            # summarise, fit and plot
            det_models: dict[str, dict] = {}
            det_points: dict[str, dict] = {}
            det_reso: dict[str, dict] = {}
            for st, per_energy in samples.items():
                aoe_samples = {
                    e: (np.concatenate(v) if len(v) > 0 else np.array([], dtype=float))
                    for e, v in per_energy.items()
                }
                model, points, stats, resolution = extract_energy_dependence(
                    aoe_samples, f"{det} in {runid} ({st})", log
                )
                if model is not None:
                    det_models[st] = model
                det_points[st] = points
                if resolution is not None:
                    det_reso[st] = resolution

                _plot_distributions(
                    pdf,
                    aoe_samples,
                    stats,
                    f"{det} in {runid}: raw A/E of electrons "
                    f"({_SIM_TYPE_STYLE[st]['label']})",
                )

            _plot_modes(
                pdf,
                det_points,
                det_models,
                f"{det} in {runid}: A/E mean energy dependence",
            )

            if det_models:
                models[det] = det_models
            stats_out[det] = {"electron_gun": det_points, "mc_resolution": det_reso}

    log.info("saving the A/E mean models")
    pars_file = Path(args.pars_file)
    pars_file.parent.mkdir(parents=True, exist_ok=True)
    dbetto.utils.write_dict(models, str(pars_file))

    stats_file = Path(args.stats_file)
    stats_file.parent.mkdir(parents=True, exist_ok=True)
    dbetto.utils.write_dict(stats_out, str(stats_file))


if __name__ == "__main__":
    main()
