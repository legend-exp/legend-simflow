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

from lgdo import Struct, lh5
from reboost import units

from legendsimflow import psl

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


ALIGNMENT_IDX = 3000  # Index to align current waveforms to Amax
NSAMPLES_OUTPUT_CURRENT_WFS = (
    5001  # Final length of the realistic current waveforms in the map
)
DT_DATA = 16.0  # Time step of the original data waveforms in ns (used to scale the derivative)
MW_PARS = {
    "length": 48,
    "num_mw": 3,
    "mw_type": 0,
}  # Parameters for the moving window average step


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--detector", required=True, help="Detector name (LH5 group)")
    parser.add_argument(
        "--sigma-conv",
        type=float,
        required=True,
        help="Sigma of the gaussian component of the convolution kernel in ns",
    )
    parser.add_argument(
        "--tau-conv",
        type=float,
        required=True,
        help="Tau of the exponential component of the convolution kernel in ns",
    )
    parser.add_argument("--input-file", required=True)
    parser.add_argument("--output-file", required=True)
    args = parser.parse_args()

    # 1. Load data
    ideal_map_obj = lh5.read(args.detector, args.input_file)
    dt = ideal_map_obj["dt"].value * units.units_convfact(ideal_map_obj["dt"], "ns")

    # 2. Setup Physics Kernel (mu=0, sigma=sigma-conv ns, tau=tau-conv ns)
    rf_kernel = psl.build_electronics_response_kernel(
        dt,
        mu_bandwidth=0,
        sigma_bandwidth=args.sigma_conv,
        tau_rc=args.tau_conv,
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
