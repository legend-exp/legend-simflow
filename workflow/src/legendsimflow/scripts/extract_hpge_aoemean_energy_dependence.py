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
the ``hit`` tier, without applying any energy-dependence correction. Only the
events that stopped the whole electron are kept: when a bremsstrahlung photon
escapes the detector the deposit shrinks while the current noise does not, and
A/E stops measuring the band. Electrons are single-site events, so what is left
is a peak with a low-side tail, fitted with the same shape as on the data side;
its position is the mean of the Gaussian. The positions are fitted as a
function of the energy with the same linear model and fitting procedure
(:class:`pygama.pargen.AoE_cal.Pol1`), which is the A/E mean model consumed by
the ``hit`` tier (see :func:`legendsimflow.hpge_pars.build_aoe_mean_func_dict`).

The dead layer is not modelled here: the A/E of an event is a ratio in which the
charge-collection efficiency largely cancels, and the correction is only needed
to centre the band.
"""

import argparse
import logging
from collections.abc import Sequence
from pathlib import Path

import awkward as ak
import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
import lh5
import numpy as np
import pyg4ometry
import pygeomtools
from dbetto import AttrsDict
from dbetto.utils import load_dict
from iminuit import Minuit, cost
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from numpy.typing import NDArray
from pygama.math.distributions import exgauss, gaussian
from pygama.pargen.AoE_cal import Pol1, aoe_peak
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, patterns, utils
from legendsimflow import reboost as reboost_utils
from legendsimflow.metadata import get_tier_settings
from legendsimflow.plot import decorate, plot_aoe_distributions, plot_aoe_vs_energy
from legendsimflow.scripts import log_script_invocation

log = logging.getLogger(__name__)


def _fit_aoe_peak(
    aoe: NDArray, nbins: int
) -> tuple[dict[str, float | int], tuple[NDArray, NDArray] | None]:
    """Position and width of the A/E peak simulated at one electron energy.

    Fits the model used on LEGEND-200 data
    (``pygama.pargen.AoE_cal.aoe_peak``): a Gaussian, the single-site peak,
    plus an exponentially modified Gaussian sharing its mean and width, the
    low-side tail left by the electrons whose charge arrives over a longer
    drift-time spread. The position is the mean of the Gaussian, which is what
    :meth:`pygama.pargen.AoE_cal.CalAoE.energy_correction` takes as the band
    centre on the data side.

    The likelihood is binned, unlike
    :func:`pygama.pargen.AoE_cal.unbinned_aoe_fit`: with thousands of events
    per fit the information lost to binning is negligible (the uncertainties
    agree to the digit) while the fit is two orders of magnitude cheaper. The
    seeding follows the unbinned routine, fitting the core and then the tail
    before the full model.

    Returns
    -------
    stats
        Keys ``n_events``, ``mu``, ``mu_err``, ``sigma``, ``chi2_ndf``. All but
        the first are NaN if the fit does not converge.
    curve
        The bin edges and the expected counts of the best fit, to draw it over
        the data. ``None`` if the fit did not converge.
    """
    aoe = np.asarray(aoe, dtype=float)
    aoe = aoe[np.isfinite(aoe)]
    n = int(aoe.size)
    failed = {
        "n_events": n,
        "mu": np.nan,
        "mu_err": np.nan,
        "sigma": np.nan,
        "chi2_ndf": np.nan,
    }
    if n < 100:
        return failed, None

    # a robust first look, to place the fit window
    mu = float(np.median(aoe))
    sigma = 1.4826 * float(np.median(np.abs(aoe - mu)))
    if sigma <= 0:
        return failed, None

    # the same window as the data side: the tail down to 15 sigma, 5 above
    lo, hi = mu - 15 * sigma, mu + 5 * sigma
    counts, edges = np.histogram(aoe, bins=nbins, range=(lo, hi))

    # widths are bounded well away from zero, in units of the data: pygama's
    # exgauss divides by the tail scale, and Minuit walks onto its bounds when
    # a tail is poorly constrained
    floor = 1e-3 * sigma

    # the Gaussian core alone, for the mean and the width
    core = (edges[:-1] >= mu - 0.5 * sigma) & (edges[1:] <= mu + 3 * sigma)
    m = Minuit(
        cost.ExtendedBinnedNLL(
            counts[core],
            np.append(edges[:-1][core], edges[1:][core][-1]),
            gaussian.cdf_ext,
        ),
        area=counts[core].sum(),
        mu=mu,
        sigma=sigma,
    )
    m.limits["sigma"] = (floor, None)
    m.migrad()
    mu, sigma = m.values["mu"], m.values["sigma"]

    # then the tail alone, below 5 sigma, for its scale and area
    tail = edges[1:] <= mu - 5 * sigma
    m = Minuit(
        cost.ExtendedBinnedNLL(
            counts[tail],
            np.append(edges[:-1][tail], edges[1:][tail][-1]),
            exgauss.cdf_ext,
        ),
        area=max(counts[tail].sum(), 1.0),
        mu=mu,
        sigma=sigma,
        tau=10 * sigma,
    )
    m.fixed["mu"] = m.fixed["sigma"] = True
    m.limits["area"] = (0, None)
    m.limits["tau"] = (floor, None)
    m.simplex().migrad()
    tau, n_bkg = m.values["tau"], m.values["area"]

    def cdf(xe, n_sig, mu, sigma, n_bkg, tau):
        return aoe_peak.cdf_ext(xe, lo, hi, n_sig, mu, sigma, n_bkg, tau)

    m = Minuit(
        cost.ExtendedBinnedNLL(counts, edges, cdf),
        n_sig=max(counts.sum() - n_bkg, 1.0),
        mu=mu,
        sigma=sigma,
        n_bkg=n_bkg,
        tau=tau,
    )
    for par in ("n_sig", "n_bkg"):
        m.limits[par] = (0, None)
    m.limits["sigma"] = (floor, None)
    m.limits["tau"] = (floor, None)
    m.migrad()
    m.hesse()

    if not m.valid or not np.isfinite(m.errors["mu"]) or m.errors["mu"] <= 0:
        return failed, None

    filled = counts > 0
    expected = np.diff(cdf(edges, *m.values))
    chi2 = float(np.sum((counts[filled] - expected[filled]) ** 2 / counts[filled]))

    return {
        "n_events": n,
        "mu": float(m.values["mu"]),
        "mu_err": float(m.errors["mu"]),
        "sigma": float(m.values["sigma"]),
        "chi2_ndf": chi2 / (int(filled.sum()) - m.nfit),
    }, (edges, expected)


def fit_aoe_peak(
    aoe: NDArray, nbins: int = 400
) -> tuple[dict[str, float | int], tuple[NDArray, NDArray] | None]:
    """:func:`_fit_aoe_peak`, with a failed fit reported instead of raised.

    A detector whose A/E at one energy defeats the fit only loses that point:
    the run it belongs to still gets its model from the other energies.
    """
    try:
        return _fit_aoe_peak(aoe, nbins)
    except Exception:
        log.warning(
            "the A/E peak fit failed on a sample of %d", len(aoe), exc_info=True
        )
        n = int(np.sum(np.isfinite(aoe)))
        return {
            "n_events": n,
            "mu": np.nan,
            "mu_err": np.nan,
            "sigma": np.nan,
            "chi2_ndf": np.nan,
        }, None


def fit_aoe_vs_energy(
    energies: Sequence[float], positions: Sequence[float], errors: Sequence[float]
) -> dict | None:
    """Fit the A/E energy dependence with the linear model used on data.

    Mirrors the mean fit of
    :meth:`pygama.pargen.AoE_cal.CalAoE.energy_correction`: the
    :class:`pygama.pargen.AoE_cal.Pol1` model (``x*a+b``) and its initial guess,
    a least-squares cost with the soft-L1 loss, and the simplex, migrad and
    hesse Minuit steps. Energies whose peak fit did not converge, and so carry
    a non-finite or vanishing uncertainty, are discarded.

    Returns
    -------
    dict | None
        The model, with the ``expression`` of the fitted function and its
        best-fit ``pars`` and ``errs``. ``None`` when fewer than two energies
        have a valid peak fit, or when Minuit does not converge.
    """
    x = np.asarray(energies, dtype=float)
    y = np.asarray(positions, dtype=float)
    ey = np.asarray(errors, dtype=float)

    # a point with no uncertainty (a distribution with no width) says nothing
    # about how well the model fits, and would dominate the least squares
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(ey) & (ey > 0)
    x, y, ey = x[valid], y[valid], ey[valid]

    if x.size < 2:
        return None

    c = cost.LeastSquares(x, y, ey, Pol1.func)
    c.loss = "soft_l1"
    m = Minuit(c, *Pol1.guess(x, y, ey))
    m.simplex()
    m.migrad()
    m.hesse()

    if not m.valid:
        return None

    return {
        "expression": Pol1.string_func("x"),
        "pars": {"a": float(m.values["a"]), "b": float(m.values["b"])},
        "errs": {"a": float(m.errors["a"]), "b": float(m.errors["b"])},
    }


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
            "electron-gun stp tier files; the electron energy is read from the "
            "file name"
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
    runid = args.runid

    log = ldfs.utils.build_log(config.metadata.simprod.config.logging, args.log_file)
    log_script_invocation(log, "extract-hpge-aoemeanmod", parser, args)

    tier_hit_settings = get_tier_settings(config, "hit")
    simulate_psd = tier_hit_settings.get("simulate_psd", True)
    simulate_psd_with_psl = tier_hit_settings.get("simulate_psd_with_psl", False)

    # the modelable detectors of this run
    is_modelable = load_dict(nersc.dvs_ro(config, args.is_modelable_file))
    dets = sorted(d for d, flag in is_modelable.get(runid, {}).items() if flag)
    log.info("modelable HPGe detectors in %s: %s", runid, ", ".join(dets))

    log.debug("loading the geometry")
    registry = pyg4ometry.gdml.Reader(
        nersc.dvs_ro(config, args.geom_file)
    ).getRegistry()
    germanium = pygeomtools.get_all_sensvols(registry, "germanium")

    # group the input files by the energy encoded in their name
    files_by_energy: dict[int, list[Path]] = {}
    for f in args.electron_stp_files:
        energy = patterns.electron_gun_energy_from_path(f)
        files_by_energy.setdefault(energy, []).append(nersc.dvs_ro(config, Path(f)))
    files_by_energy = {e: sorted(files_by_energy[e]) for e in sorted(files_by_energy)}
    log.info(
        "electron energies found: %s keV", ", ".join(str(e) for e in files_by_energy)
    )

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

            if det not in germanium:
                log.warning(
                    "%s is not a germanium volume of the geometry, skipping", det
                )
                continue

            # the detector origin is where the geometry places its volume
            det_loc = registry.physicalVolumeDict[det].position

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

            # collect the raw A/E values per PSD method and electron energy.
            # one detector at one energy is a few MB at most, so the steps are
            # read in one go
            samples: dict[str, dict[int, NDArray]] = {
                method: dict.fromkeys(files_by_energy, np.array([], dtype=float))
                for method in methods
            }
            for energy, files in files_by_energy.items():
                # a job might have recorded no hit at all in this detector
                with_det = [str(f) for f in files if lh5.ls(f, f"stp/{det}")]
                if not with_det:
                    log.warning("no %s hits found in the %d keV files", det, energy)
                    continue

                steps = lh5.read(f"stp/{det}", with_det).view_as("ak", with_units=True)
                energy_dep = ak.sum(steps.edep, axis=-1)

                # keep only the events that stopped the whole electron in this
                # detector. When a bremsstrahlung photon escapes, the deposit is
                # smaller while the current noise, which is an absolute
                # amplitude, is not: A/E is then scattered over orders of
                # magnitude and no longer measures the band. The data side has
                # no such events either, its Compton bands fix the energy
                full_energy = np.asarray(energy_dep > 0.999 * energy)

                psd, psd_detailed = reboost_utils.compute_hpge_psd_observables(
                    steps,
                    steps.edep,
                    energy_dep,
                    det_loc,
                    psd_inputs,
                    # no correction and no classifier: only the raw A/E is of
                    # interest here
                    aoe_res=1.0,
                    aoe_mean=1.0,
                    aoe_mean_psl=1.0,
                )
                for method, values in (
                    ("single_template", psd.aoe),
                    ("psl", psd_detailed.aoe),
                ):
                    if method in samples:
                        samples[method][energy] = np.asarray(values, dtype=float)[
                            full_energy
                        ]

            # summarise, fit and plot
            det_models: dict[str, dict] = {}
            det_points: dict[str, dict] = {}
            for method, aoe in samples.items():
                energies = sorted(aoe)
                fits = {energy: fit_aoe_peak(aoe[energy]) for energy in energies}
                stats = {energy: f[0] for energy, f in fits.items()}
                curves = {energy: f[1] for energy, f in fits.items()}
                det_points[method] = {
                    "energy_in_keV": energies,
                    **{
                        key: [stats[energy][key] for energy in energies]
                        for key in ("n_events", "mu", "mu_err", "sigma", "chi2_ndf")
                    },
                }

                for energy in energies:
                    log.info(
                        "%s in %s (%s) at %d keV: n = %d, mu = %.5f +- %.5f, "
                        "sigma = %.2e, chi2/ndf = %.1f",
                        det,
                        runid,
                        method,
                        energy,
                        *(
                            stats[energy][key]
                            for key in ("n_events", "mu", "mu_err", "sigma", "chi2_ndf")
                        ),
                    )

                model = fit_aoe_vs_energy(
                    energies,
                    det_points[method]["mu"],
                    det_points[method]["mu_err"],
                )
                if model is None:
                    log.warning(
                        "%s in %s (%s): no A/E energy-dependence model, the fit did "
                        "not converge or too few energies have a valid peak fit",
                        det,
                        runid,
                        method,
                    )
                else:
                    det_models[method] = model
                    log.info(
                        "%s in %s (%s): a = %.3e +- %.1e 1/keV, b = %.5f +- %.1e",
                        det,
                        runid,
                        method,
                        model["pars"]["a"],
                        model["errs"]["a"],
                        model["pars"]["b"],
                        model["errs"]["b"],
                    )

                fig = plot_aoe_distributions(
                    aoe,
                    stats,
                    curves,
                    f"{det} in {runid}: raw A/E of electrons ({method})",
                )
                decorate(fig)
                pdf.savefig(fig)
                plt.close(fig)

            fig = plot_aoe_vs_energy(
                det_points, det_models, f"{det} in {runid}: A/E mean energy dependence"
            )
            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

            if det_models:
                models[det] = det_models
            stats_out[det] = det_points

    log.info("saving the A/E mean models")
    for path, content in ((args.pars_file, models), (args.stats_file, stats_out)):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        dbetto.utils.write_dict(content, str(path))


if __name__ == "__main__":
    main()
