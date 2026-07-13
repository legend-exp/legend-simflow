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

import argparse
import contextlib
import logging
import warnings
from collections.abc import Iterator
from pathlib import Path

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
import numpy as np
from iminuit import Minuit
from lgdo import lh5
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pygama.pargen.AoE_cal import CalAoE
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import utils
from legendsimflow.plot import decorate
from legendsimflow.scripts import log_script_invocation

DEFAULT_SETTINGS = {"simid_regex": "*Pb212*"}


def get_aoe_reso(m: Minuit, x: float = 2040) -> tuple[float, float]:
    """Get the A/E resolution and its uncertainty at a given energy x from a Minuit fit object m.

    Based on simple propagation of errors.

    Parameters
    ----------
    m
        The Minuit fit object containing the fit results and covariance matrix.
    x
        The energy at which to evaluate the A/E resolution. Default is 2040 keV.

    Returns
    -------
    tuple[float, float]
        A tuple containing the A/E resolution and its uncertainty at energy x, both multiplied by 100 to convert to percentage.

    """
    if not getattr(m, "valid", False):
        msg = "Minuit fit is not valid (m.valid is False)."
        raise RuntimeError(msg)
    if getattr(m, "covariance", None) is None:
        msg = "Covariance matrix is missing. Run m.hesse() after migrad()."
        raise RuntimeError(msg)

    a = m.values["a"]
    b = m.values["b"]
    c = m.values["c"]

    Caa = m.covariance["a", "a"]
    Cbb = m.covariance["b", "b"]
    Ccc = m.covariance["c", "c"]
    Cab = m.covariance["a", "b"]
    Cac = m.covariance["a", "c"]
    Cbc = m.covariance["b", "c"]

    x_arr = np.asarray(x, dtype=float)
    u = x_arr + 1e-99
    r = b / u

    # Real-domain requirements for this analytic derivative form:
    #   t = r**c and ln(r) require r > 0
    if np.any(r <= 0):
        msg = "Need b/(x+1e-99) > 0 for real-valued power/log derivatives."
        raise ValueError(msg)

    t = r**c
    g = a + t
    if np.any(g <= 0):
        msg = "Need a + (b/(x+1e-99))**c > 0 for sqrt."
        raise ValueError(msg)

    y = np.sqrt(g)

    # Correct derivatives
    dfa = 1.0 / (2.0 * y)
    dfb = (c * t) / (2.0 * y * b)
    dfc = (t * np.log(r)) / (2.0 * y)

    var = (
        Caa * dfa**2
        + Cbb * dfb**2
        + Ccc * dfc**2
        + 2.0 * Cab * dfa * dfb
        + 2.0 * Cac * dfa * dfc
        + 2.0 * Cbc * dfb * dfc
    )

    sigma = np.sqrt(np.clip(var, 0.0, None))

    if np.ndim(x_arr) == 0:
        return 100 * float(y), 100 * float(sigma)
    return 100 * y, 100 * sigma


def get_mean_err(m: Minuit, x=2040):
    """Get the A/E mean and its uncertainty at a given energy x from a Minuit fit object m.

    Based on simple propagation of errors.

    Parameters
    ----------
    m
        The Minuit fit object containing the fit results and covariance matrix.
    x
        The energy at which to evaluate the A/E mean. Default is 2040 keV.

    Returns
    -------
    tuple[float, float]
        A tuple containing the A/E mean and its uncertainty at energy x, both multiplied by 100 to convert to percentage.

    """
    if not getattr(m, "valid", False):
        msg = "Minuit fit is not valid (m.valid is False)."
        raise RuntimeError(msg)
    if getattr(m, "covariance", None) is None:
        msg = "Covariance matrix is missing. Run m.hesse() after migrad()."
        raise RuntimeError(msg)

    a = m.values["a"]
    b = m.values["b"]

    # covariance entries
    var_a = m.covariance["a", "a"]  # sigma_a^2
    var_b = m.covariance["b", "b"]  # sigma_b^2
    cov_ab = m.covariance["a", "b"]  # covariance(a, b)

    # evaluate at particular x0
    y0 = b + a * x

    # propagated variance: var(y0) = var_b + x0^2*var_a + 2*x0*cov_ab
    var_y0 = var_b + x**2 * var_a + 2 * x * cov_ab

    return 100 * y0, 100 * np.sqrt(var_y0)


@contextlib.contextmanager
def suppress_fit_warnings() -> Iterator[None]:
    """Silence the harmless numpy RuntimeWarnings emitted during the A/E fits.

    The compton-band and energy-dependence fits routinely explore parameter
    regions that emit RuntimeWarnings (e.g. sqrt of a negative value). These are
    suppressed locally so they are not promoted to errors (e.g. by a
    ``filterwarnings = error`` test configuration), which would otherwise abort
    the fit inside pygama.

    Yields
    ------
    None

    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        yield


def get_mc_reso_entry(cal: CalAoE, log: logging.Logger) -> dict[str, float | None]:
    """A/E resolution entry for a fitted :class:`CalAoE`.

    Returns a ``{"val": ..., "err": ...}`` dictionary with the A/E resolution at
    2040 keV (in percent) and its uncertainty. If the sigma-vs-energy fit did not
    converge to a valid result, both values are ``None``: the mean energy
    correction is the primary product of this script, while the resolution is a
    quality-assurance quantity, so a non-convergent resolution fit must not abort
    the whole computation.

    Parameters
    ----------
    cal
        A :class:`CalAoE` instance on which ``energy_correction()`` has been run.
    log
        Logger used to report a non-convergent fit.

    Returns
    -------
    dict[str, float | None]
        The resolution value and uncertainty, or ``None`` values on failure.

    """
    sigma_fit = getattr(cal, "SigmaFit_obj", None)
    if sigma_fit is not None:
        try:
            with suppress_fit_warnings():
                val, err = get_aoe_reso(sigma_fit)
            return {"val": float(val), "err": float(err)}
        except (RuntimeError, ValueError) as e:
            log.warning("A/E resolution fit did not converge: %s", e)
    return {"val": None, "err": None}


@snakemake_compatible(
    mapping={
        "runid": "wildcards.runid",
        "hpge_detector": "wildcards.hpge_detector",
        "pars_file": "output.pars_file",
        "plot_file": "output.plot_file",
        "settings": "input.settings",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main():
    parser = argparse.ArgumentParser(
        description="Compute A/E energy-dependence correction parameters for an HPGe detector."
    )

    parser.add_argument(
        "--hpge-detector",
        required=True,
        dest="hpge_detector",
        help="HPGe detector name",
    )
    parser.add_argument(
        "--pars-file",
        required=True,
        help="output YAML file for A/E energy-correction parameters (and validation results)",
    )
    parser.add_argument(
        "--settings",
        type=str,
        required=False,
        default=None,
        help="Path to YAML file with settings for the fitting ",
    )
    parser.add_argument(
        "--plot-file",
        required=True,
        help="output PDF plot file",
    )
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
    log_script_invocation(log, "get-aoe-energy-correction", parser, args)

    det = args.hpge_detector

    settings = (
        dbetto.AttrsDict(dbetto.utils.load_dict(args.settings))
        if args.settings is not None
        else dbetto.AttrsDict(DEFAULT_SETTINGS)
    )

    # only if we find ssc data
    fit_data = False

    if (l200data := Path(config.paths.l200data)).exists():
        data_files = sorted((l200data / "generated/tier/pht/ssc/p16").glob("*/*.lh5"))
    else:
        data_files = []

    fit_data = len(data_files) > 0

    if fit_data:
        log.info("... starting reading data")

        l200data = Path(config.paths.l200data)

        data = lh5.read(
            f"hit/{det}",
            data_files,
            field_mask=["AoE_Corrected", "cuspEmax_ctc_cal"],
        ).view_as("pd")

        is_pulser = (
            lh5.read(
                "ch1027201/dsp/trapTmax",
                [Path(str(p).replace("/pht/", "/psp/")) for p in data_files],
            ).view_as("pd")
            > 200
        )

        data = data[(data.cuspEmax_ctc_cal > 600) & (~is_pulser)]

        log.info("... read data")

        aoe_cal = CalAoE(compt_bands_width=50)
        with suppress_fit_warnings():
            aoe_cal.energy_correction(
                data,
                aoe_param="AoE_Corrected",
                corrected_param="AoE_reCorrected",
                display=0,
            )
        log.info("... fitted data")

        data_mu, data_mu_err = get_mean_err(aoe_cal.mean_fit_obj)
        data_sig, data_sig_err = get_aoe_reso(aoe_cal.SigmaFit_obj)

    log.info(" ... load MC")

    path_mc = config.paths.tier_hit

    files_mc = list(
        Path(path_mc).glob(f"generated/tier/hit/{settings.simid_regex}/*.lh5")
    )

    mc = lh5.read(
        f"hit/{det}",
        files_mc,
        field_mask=[
            "psd/single_temp/aoe_raw",
            "psd/pulse_lib/aoe_raw",
            "energy",
        ],
    ).view_as("pd")

    log.info("... fitting mc")
    aoe_cal_mc = CalAoE(compt_bands_width=50, cal_energy_param="energy")
    with suppress_fit_warnings():
        aoe_cal_mc.energy_correction(
            mc,
            aoe_param="psd_single_temp_aoe_raw",
            corrected_param="AoE_reCorrected",
            display=0,
        )

    aoe_cal_mc_d = CalAoE(compt_bands_width=50, cal_energy_param="energy")
    with suppress_fit_warnings():
        aoe_cal_mc_d.energy_correction(
            mc,
            aoe_param="psd_pulse_lib_aoe_raw",
            corrected_param="AoE_reCorrected",
            display=0,
        )

    log.info("... saving pars")
    mfs = aoe_cal_mc.energy_corr_res_dict["mean_fits"]
    mfs_d = aoe_cal_mc_d.energy_corr_res_dict["mean_fits"]

    out = {
        "single_template": {
            "expression": mfs["expression"],
            "pars": mfs["pars"],
        },
        "psl": {"expression": mfs_d["expression"], "pars": mfs_d["pars"]},
    }
    mc_reso = {
        "single_template": get_mc_reso_entry(aoe_cal_mc, log),
        "psl": get_mc_reso_entry(aoe_cal_mc_d, log),
    }

    data_out = None
    if fit_data:
        data_out = {
            "mu": {"val": float(data_mu), "err": float(data_mu_err)},
            "sigma": {"val": float(data_sig), "err": float(data_sig_err)},
        }

    output = {"energy_corrections": out, "mc_resolution": mc_reso}
    if fit_data and data_out is not None:
        output["data_validation"] = data_out

    pars_file = Path(args.pars_file)
    pars_file.parent.mkdir(parents=True, exist_ok=True)
    dbetto.utils.write_dict(output, str(pars_file))

    plot_file = Path(args.plot_file)
    plot_file.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(str(plot_file)) as pdf:
        log.info("... plotting")

        x = aoe_cal_mc.energy_corr_fits["mean"].index
        y = aoe_cal_mc.energy_corr_fits["mean"].to_numpy()
        yerr = aoe_cal_mc.energy_corr_fits["mean_err"].to_numpy()

        fig, ax = plt.subplots(figsize=(12, 3.5))

        ax.errorbar(
            x,
            y,
            yerr=yerr,
            capsize=1,
            marker=".",
            color="#0077BB",
            label="MC (single-temp)",
            linestyle="none",
            markersize=10,
        )

        ax.plot(
            np.linspace(900, 2400, 1000),
            aoe_cal_mc.mean_func.func(
                np.linspace(900, 2400, 1000),
                **aoe_cal_mc.energy_corr_res_dict["mean_fits"]["pars"],
            ),
            color="#1A2A5B",
            label=f"Best fit (st) {mfs['pars']['a'] * 10**5:.2f} %/MeV",
        )

        x = aoe_cal_mc_d.energy_corr_fits["mean"].index
        y = aoe_cal_mc_d.energy_corr_fits["mean"].to_numpy()
        yerr = aoe_cal_mc_d.energy_corr_fits["mean_err"].to_numpy()

        ax.errorbar(
            x,
            y,
            yerr=yerr,
            capsize=1,
            marker=".",
            linestyle="none",
            color="darkorange",
            label="MC (pulse lib)",
            markersize=10,
        )

        ax.plot(
            np.linspace(900, 2400, 1000),
            aoe_cal_mc_d.mean_func.func(
                np.linspace(900, 2400, 1000),
                **aoe_cal_mc_d.energy_corr_res_dict["mean_fits"]["pars"],
            ),
            color="#CC3311",
            label=f"Best fit (psl) {mfs_d['pars']['a'] * 10**5:.2f} %/MeV",
        )

        ax.set_xlabel("energy [keV]")
        ax.set_ylabel("A/E Corrected")
        ax.legend()

        ax.set_title(f"Energy dependence in SSC MC {det}")

        decorate(fig)
        pdf.savefig()
        plt.close(fig)

        x = aoe_cal_mc.energy_corr_fits["sigma"].index
        y = 100 * aoe_cal_mc.energy_corr_fits["sigma"].to_numpy()
        yerr = 100 * aoe_cal_mc.energy_corr_fits["sigma_err"].to_numpy()

        fig, ax = plt.subplots(figsize=(12, 3.5))
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            capsize=1,
            marker=".",
            color="#0077BB",
            linestyle="none",
            markersize=10,
            label="MC",
        )

        e = np.linspace(900, 2400, 1000)
        ax.plot(
            e,
            100
            * aoe_cal_mc.sigma_func.func(
                e, **aoe_cal_mc.energy_corr_res_dict["SigmaFits"]["pars"]
            ),
            color="#1A2A5B",
            label="Best fit (MC)",
        )

        if fit_data:
            x = aoe_cal.energy_corr_fits["sigma"].index
            y = 100 * aoe_cal.energy_corr_fits["sigma"].to_numpy()
            yerr = 100 * aoe_cal.energy_corr_fits["sigma_err"].to_numpy()

            ax.errorbar(
                x,
                y,
                yerr=yerr,
                capsize=1,
                marker=".",
                linestyle="none",
                markersize=10,
                color="darkorange",
                label="Data",
            )

            ax.plot(
                e,
                100
                * aoe_cal.sigma_func.func(
                    e, **aoe_cal.energy_corr_res_dict["SigmaFits"]["pars"]
                ),
                color="#CC3311",
                label="Best fit (data)",
            )

        ax.set_xlabel("energy [keV]")
        ax.set_ylabel("A/E resolution [%]")
        ax.legend(ncol=2)
        ax.set_title(f"Energy dependence in SSC data {det}")

        ax.set_ylim(0.2, 1.5)

        decorate(fig)

        pdf.savefig()
        plt.close(fig)


if __name__ == "__main__":
    main()
