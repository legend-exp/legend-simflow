from __future__ import annotations

import argparse

import dbetto

from legendsimflow.hpge_pars import get_pulse_pars


def parse_args():
    parser = argparse.ArgumentParser(description="Example script")

    parser.add_argument(
        "-i",
        "--data-path",
        required=True,
        help="Path to the data production environment",
    )
    parser.add_argument(
        "-p",
        "--plot-path",
        required=False,
        default=None,
        help="Path to the folder to save plots.",
    )
    parser.add_argument(
        "-o", "--out-file", required=True, help="Path to the output file"
    )

    parser.add_argument(
        "-r", "--run", type=str, required=True, help="Run number formatted as pXY-rXYZ"
    )
    parser.add_argument(
        "-d", "--detectors", nargs="+", required=True, help="One or more detector names"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    period = args.run.split("-")[0]
    run = args.run.split("-")[1]

    out = {
        det: get_pulse_pars(
            det, args.data_path, period=period, run=run, plot_path=args.plot_path
        )
        for det in args.detectors
    }

    dbetto.utils.write_dict(out, args.out_file)
