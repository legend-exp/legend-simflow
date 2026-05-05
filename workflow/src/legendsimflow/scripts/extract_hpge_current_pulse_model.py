# ruff: noqa: I002

# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>,
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

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from reboost.hpge.psd import _current_pulse_model as current_pulse_model
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import hpge_pars, utils
from legendsimflow import metadata as mutils
from legendsimflow.plot import decorate
from legendsimflow.scripts import log_script_invocation


def _plot_currmod_from_metadata(
    entry: dict, hpge: str, runid: str, plot_file: str
) -> None:
    pars = entry["current_pulse_pars"]
    t = np.linspace(-500, 1500, 2000)
    y = current_pulse_model(t, **pars)
    fig, ax = plt.subplots()
    ax.plot(t, y)
    ax.set_xlim(-500, 1500)
    ax.set_xlabel("time [ns]")
    ax.set_ylabel("current [a.u.]")
    ax.set_title(f"{hpge} in {runid}: current pulse model (from metadata)")
    with PdfPages(plot_file) as pdf:
        pdf.savefig(fig)
    plt.close(fig)


@snakemake_compatible(
    mapping={
        "runid": "wildcards.runid",
        "hpge_detector": "wildcards.hpge_detector",
        "pars_file": "output.pars_file",
        "plot_file": "output.plot_file",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract the HPGe current pulse model for a LEGEND run."
    )
    parser.add_argument(
        "--runid",
        required=True,
        help="LEGEND run identifier (e.g. l200-p03-r000-phy)",
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
        help="output YAML file for the current pulse model parameters",
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
    log_script_invocation(log, "extract-hpge-currmod", parser, args)

    runid = args.runid
    hpge = args.hpge_detector
    pars_file = args.pars_file
    plot_file = args.plot_file

    # use just a single waveform
    single_waveform = True

    # check for metadata-driven defaults first; if present, bypass the
    # waveform fitting entirely (enables LEGEND-1000 simulations without
    # l200data)
    raw_currmod = mutils.simpars(
        metadata, "geds.currmod", runid, config.experiment, default=None
    )
    currmod_default = (
        raw_currmod.get("default", None) if raw_currmod is not None else None
    )

    if currmod_default is not None:
        log.info("... using currmod metadata defaults for %s in %s", hpge, runid)
        entry = raw_currmod.get(hpge, currmod_default)
        dbetto.utils.write_dict(entry.to_dict(), pars_file)
        _plot_currmod_from_metadata(entry.to_dict(), hpge, runid, plot_file)
        return

    l200data = getattr(config.paths, "l200data", None)
    if l200data is None:
        msg = (
            "l200data is not configured and no currmod metadata with a 'default' "
            "key was found — cannot extract current pulse model"
        )
        raise RuntimeError(msg)

    hit_tier_name = utils.get_hit_tier_name(l200data)

    msg = f"... determined hit tier name is {hit_tier_name}"
    log.info(msg)
    log.info("... looking up the fit inputs")

    msg = f"... using a single waveform? {single_waveform}"
    log.info(msg)

    raw_wf_pairs, dsp_cfg_file, all_dts, selected_dts = (
        hpge_pars.lookup_currmod_fit_inputs(
            l200data,
            metadata,
            runid,
            hpge,
            hit_tier_name,
            max_waveforms=1 if single_waveform else 100,
        )
    )

    lh5_group = mutils._get_lh5_table(
        metadata,
        raw_wf_pairs[0][0],
        hpge,
        "raw",
        runid,
    )

    log.info("... fetching %d current pulse(s)", len(raw_wf_pairs))
    times_list, current_list = hpge_pars.get_current_pulses(
        raw_wf_pairs, lh5_group, str(dsp_cfg_file)
    )

    log.info("... fitting the current pulse(s) to extract the model")
    popt, x, y = hpge_pars.fit_currmod(times_list, current_list)

    # now plot
    log.info("... plotting the fit result")

    with PdfPages(plot_file) as pdf:
        log.info("... plotting drift-time selection")

        fig, _ = hpge_pars.plot_dt_selection(all_dts, selected_dts)
        fig.suptitle(f"{hpge} in {runid}: drift-time selection for current-pulse fit")
        decorate(fig)
        pdf.savefig()
        plt.close(fig)
        del all_dts, selected_dts

        fig, ax = hpge_pars.plot_currmod_fit_result(times_list, current_list, x, y)

        ax.set_xlim(-1200, 1200)
        fig.suptitle(f"{hpge} in {runid}: current waveform fit result")
        decorate(fig)
        pdf.savefig()

        # also save with a narrower range
        ax.set_xlim(-400, 300)
        pdf.savefig()
        plt.close(fig)
        del times_list, current_list

        log.info("... adding the mean aoe")

        popt_dict = utils._curve_fit_popt_to_dict(popt)
        mean_aoe = hpge_pars.estimate_mean_aoe(popt)

        log.info("... estimating effect of noise")

        files = hpge_pars.lookup_file_paths(
            l200data, runid, hit_tier_name=hit_tier_name
        )

        temp = current_pulse_model(np.linspace(-500, 1000, 1501), *popt)

        noise_sample, a_max = hpge_pars.get_noise_maxima_and_sample(
            files.raw,
            files.hit,
            lh5_group,
            str(dsp_cfg_file),
            "curr_av",
            temp,
            norm=mean_aoe * 2000,
            sample_size=100,
            maximum_number=100000,
        )

        log.info("... plot noise waveforms")
        fig, ax = hpge_pars.plot_noise_waveforms(
            noise_sample, temp, norm=mean_aoe * 2000
        )
        fig.suptitle(f"{hpge} in {runid}: noise waveforms")
        decorate(fig)
        pdf.savefig()
        plt.close(fig)
        del noise_sample, temp

        # now do the plot
        fit_result = hpge_pars.fit_noise_gauss(a_max, bins=1000)
        fig, ax = hpge_pars.plot_gauss_fit(
            a_max, fit_result, nominal_val=mean_aoe * 2000
        )
        fig.suptitle(f"{hpge} in {runid}: noise waveforms fit result")
        decorate(fig)
        pdf.savefig()
        plt.close(fig)

        log.info("... saving outputs")
        dbetto.utils.write_dict(
            {
                "current_pulse_pars": popt_dict,
                "mean_aoe": mean_aoe,
                "current_reso": fit_result.values["sigma"],
            },
            pars_file,
        )


if __name__ == "__main__":
    main()
