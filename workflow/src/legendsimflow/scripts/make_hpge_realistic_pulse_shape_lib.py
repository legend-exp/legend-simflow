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
from lgdo import Struct
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from reboost import units
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import psl, utils
from legendsimflow.plot import decorate
from legendsimflow.scripts import log_script_invocation

ALIGNMENT_IDX = 3000  # Index to align current waveforms to Amax
NSAMPLES_OUTPUT_CURRENT_WFS = (
    5001  # Final length of the realistic current waveforms in the map
)
DT_DATA = (
    psl.DT_DATA
)  # Time step of the original data waveforms in ns (used to scale the derivative)
MW_PARS = psl.MW_PARS  # Parameters for the moving window average step


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
        "--plot-file", required=False, default = None, help="Path to save validation plots"
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

    log_file = args.log_file

    log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "realistic-psl", parser, args)

    electronics_model = dbetto.utils.load_dict(args.electronics_model_file)
    if args.detector not in electronics_model:
        msg = f"Detector {args.detector} not found in '{args.electronics_model_file}'"
        raise KeyError(msg)
    try:
        detector_model = electronics_model[args.detector]
        sigma_conv = detector_model["sigma"]
        tau_conv = detector_model["tau"]
    except KeyError as e:
        missing_key = str(e)
        msg = (
            f"missing key {missing_key} in electronics-model parameters for detector "
            f"{args.detector} in {args.electronics_model_file}"
        )
        raise KeyError(msg) from e

    # 1. Load data
    ideal_map_obj = lh5.read(args.detector, args.input_file)
    dt = ideal_map_obj["dt"].value * units.units_convfact(ideal_map_obj["dt"], "ns")

    # 2. Setup Physics Kernel (mu=0, sigma=sigma ns, tau=tau ns)
    rf_kernel = psl.build_electronics_response_kernel(
        dt,
        mu_bandwidth=0,
        sigma_bandwidth=sigma_conv,
        tau_rc=tau_conv,
    )

    # 3. Process
    realistic_dict = psl.make_realistic_pulse_shape_lib(
        ideal_map_obj,
        rf_kernel,
        ALIGNMENT_IDX,
        NSAMPLES_OUTPUT_CURRENT_WFS,
        mw_pars=MW_PARS,
        dt_data=DT_DATA,
    )

    # 4. Write output with units
    out_struct = Struct(realistic_dict)
    lh5.write(
        obj=out_struct, name=args.detector, lh5_file=args.output_file, wo_mode="of"
    )

    log.info("Realistic library created successfully: %s", args.output_file)

    # 5. Validation plots
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

        log.info("validation plots saved to %s", plot_file)


if __name__ == "__main__":
    main()
