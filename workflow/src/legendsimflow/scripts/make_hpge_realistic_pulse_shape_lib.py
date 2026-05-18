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


from __future__ import annotations

import argparse
import logging

import dbetto
import lh5
from lgdo import Struct
from reboost import units
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import psl

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


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
        "output_file": "output[0]",
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
    parser.add_argument("--output-file", required=True)
    args = parser.parse_args()

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

    logger.info("Realistic library created successfully: %s", args.output_file)


if __name__ == "__main__":
    main()
