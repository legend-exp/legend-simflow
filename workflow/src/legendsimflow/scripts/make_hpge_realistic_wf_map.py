from __future__ import annotations

import argparse
import logging

from lgdo import lh5

from legendsimflow import psl

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


ALIGNMENT_IDX = 3000  # Index to align current waveforms to Amax
NSAMPLES_OUTPUT_CURRENT_WFS = (
    4001  # Final length of the realistic current waveforms in the map
)


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
    ideal_map_obj = lh5.read(f"{args.detector}/", args.input_file)
    dt = ideal_map_obj["dt"].value

    # 2. Setup Physics Kernel (mu=0, sigma=sigma-conv ns, tau=tau-conv ns)
    rf_kernel = psl.create_response_kernel(
        dt, mu_bandwidth=0, sigma_bandwidth=args.sigma_conv, tau_rc=args.tau_conv
    )

    # 3. Process
    realistic_dict = psl.generate_realistic_map(
        ideal_map_obj, rf_kernel, ALIGNMENT_IDX, NSAMPLES_OUTPUT_CURRENT_WFS
    )

    # 4. Write output with units
    psl.write_realistic_struct_to_lh5(
        realistic_dict=realistic_dict,
        original_struct=ideal_map_obj,
        output_filename=args.output_file,
        group_name=args.detector,
    )
    logger.info("Realistic library created successfully: %s", args.output_file)


if __name__ == "__main__":
    main()
