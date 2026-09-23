# ruff: noqa: I002

# Copyright (C) 2026 Giovanna Saleh <giovanna.saleh@phd.unipd.it>,
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
from pathlib import Path

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
import numpy as np
import pyg4ometry
import pygeomhpges
from lgdo import Array, Struct
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from reboost import units
from reboost.hpge import make_hpge_pulse_shape_library, plot_psl_aoe_maps
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, psl, utils
from legendsimflow.plot import decorate
from legendsimflow.scripts import log_script_invocation

ALIGNMENT_IDX = 1000  # Index to align current waveforms to Amax
NSAMPLES_OUTPUT_CURRENT_WFS = (
    4001  # Final length of the realistic current waveforms in the map
)


@snakemake_compatible(
    mapping={
        "detector": "wildcards.hpge_detector",
        "electronics_model_file": "input.electronics_model",
        "input_file": "input.ideal_psl",
        "output_file": "output.psl_file",
        "plot_file": "output.plot_file",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--detector", required=True, help="Detector name (LH5 group)")
    parser.add_argument(
        "--electronics-model-file",
        required=True,
        help=(
            "YAML file with electronics-model parameters keyed by detector; "
            "used to look up sigma/tau for the selected detector"
        ),
    )
    parser.add_argument("--input-file", required=True)
    parser.add_argument("--output-file", required=True, help="Path to output LH5 file")
    parser.add_argument(
        "--plot-file",
        required=False,
        default=None,
        help="Path to save validation plots",
    )
    parser.add_argument(
        "--dtype",
        choices=("float32", "float64"),
        default="float32",
        help="bit depth of the pulse-shape samples (default: %(default)s)",
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

    elecmod_file = nersc.dvs_ro(config, args.electronics_model_file)
    ideal_psl_file = nersc.dvs_ro(config, args.input_file)

    log_file = args.log_file

    log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "realistic-psl", parser, args)

    electronics_model = dbetto.utils.load_dict(elecmod_file)
    if args.detector not in electronics_model:
        msg = f"Detector {args.detector} not found in '{elecmod_file}'"
        raise KeyError(msg)
    try:
        detector_model = electronics_model[args.detector]
        sigma_conv = detector_model["sigma"]
        tau_conv = detector_model["tau"]
    except KeyError as e:
        missing_key = str(e)
        msg = (
            f"missing key {missing_key} in electronics-model parameters for detector "
            f"{args.detector} in {elecmod_file}"
        )
        raise KeyError(msg) from e

    # 1. Load data
    ideal_map_obj = lh5.read(args.detector, ideal_psl_file)
    dt = ideal_map_obj["dt"].value * units.units_convfact(ideal_map_obj["dt"], "ns")

    # 2. Setup Physics Kernel (mu=0, sigma=sigma ns, tau=tau ns)
    kernel_start = -100
    rf_kernel = psl.build_electronics_response_kernel(
        dt,
        mu_bandwidth=0,
        sigma_bandwidth=sigma_conv,
        tau_rc=tau_conv,
        kernel_start=kernel_start,
    )

    # 3. Process
    realistic_dict = psl.make_realistic_pulse_shape_lib(
        ideal_map_obj,
        rf_kernel,
        ALIGNMENT_IDX,
        NSAMPLES_OUTPUT_CURRENT_WFS,
        dtype=np.dtype(args.dtype),
        kernel_t0_idx=-2 * kernel_start,
    )
    # 4. normalise the current waveforms
    h_aoe, mean_aoe = psl.get_avg_aoe(
        [realistic_dict[k] for k in realistic_dict if "waveform" in k]
    )

    for key in realistic_dict:
        if "waveform" in key:
            realistic_dict[key] = Array(realistic_dict[key].view_as("np") / mean_aoe)

    # 5. Write output with units
    out_struct = Struct(realistic_dict)
    lh5.write(
        obj=out_struct, name=args.detector, lh5_file=args.output_file, wo_mode="of"
    )

    log.info("Realistic library created successfully: %s", args.output_file)

    # 6 . Validation plots
    if args.plot_file is not None:
        plot_file = Path(args.plot_file)
        plot_file.parent.mkdir(parents=True, exist_ok=True)

        angle_keys = [k for k in ideal_map_obj if "waveform" in k]
        with PdfPages(str(plot_file)) as pdf:
            for key in sorted(angle_keys):
                angle = int(key.split("_")[1])
                for scan in ("r", "z"):
                    fig, _ = psl.plot_rz_scan(
                        ideal_map_obj,
                        angle_deg=angle,
                        detector_id=args.detector,
                        scan=scan,
                        step=10,
                        xlim=(-100, 3000),
                    )
                    decorate(fig)
                    pdf.savefig(fig)
                    plt.close(fig)

                    fig, _ = psl.plot_rz_scan(
                        realistic_dict,
                        angle_deg=angle,
                        detector_id=args.detector,
                        scan=scan,
                        step=10,
                        xlim=(-1000, 1000),
                    )
                    decorate(fig)
                    pdf.savefig(fig)
                    plt.close(fig)

            fig, ax = plt.subplots()
            h_aoe.plot(ax=ax, yerr=False)

            ax.set_xlabel("A_max [arb]")
            ax.axvline(
                mean_aoe,
                color="red",
                linestyle="--",
                label=f"mean A/E = {mean_aoe:.2f}",
            )
            ax.legend()
            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

            # A/E over the (r, z) plane, one panel per crystal axis
            reg = pyg4ometry.geant4.Registry()
            natge = pygeomhpges.materials.make_natural_germanium(registry=reg)
            hpge_profile = pygeomhpges.make_hpge(
                metadata.hardware.detectors.germanium.diodes[args.detector],
                registry=reg,
                material=natge,
                allow_cylindrical_asymmetry=False,
            )

            # the waveforms are already normalised to the mean A/E
            fig, _ = plot_psl_aoe_maps(
                {
                    int(k.split("_")[1]): make_hpge_pulse_shape_library(
                        realistic_dict, k
                    )
                    for k in realistic_dict
                    if "waveform" in k
                },
                hpge=hpge_profile,
                normalise=False,
                title=args.detector,
            )
            decorate(fig)
            pdf.savefig(fig)
            plt.close(fig)

        log.info("validation plots saved to %s", plot_file)


if __name__ == "__main__":
    main()
