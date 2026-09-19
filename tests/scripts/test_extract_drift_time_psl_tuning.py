from __future__ import annotations

import sys
from pathlib import Path

import dbetto.utils
import lh5
import numpy as np
from lgdo import Array, Scalar, Struct
from pytest import fixture
from scipy.stats import norm

from legendsimflow.scripts import extract_drift_time_psl_tuning


@fixture
def make_psl_scan(tmp_path):
    out = {}
    ideal_psl = str(tmp_path / "outputs" / "l200-p16-r008-ssc-ideal_psl.lh5")

    for sidx, slope in enumerate(np.linspace(-1.0, 1.0, 3)):
        out[f"slope_{sidx}"] = {}
        for didx, depv in enumerate(np.linspace(500, 800, 3)):
            t = np.arange(5000)
            wfs = []
            # Gaussian PDF
            for mu in np.linspace(0, 1990, 200):
                wfs.append(np.cumsum(norm.pdf(t, loc=mu, scale=100)))
                wfs[-1] /= wfs[-1][-1]  # normalize to 1 at the end of the waveform

            wfs = np.array(wfs).reshape(20, 10, 5000)
            drift_time = np.linspace(0, 2000, 200).reshape(20, 10)

            r = np.linspace(0, 50, 20)
            z = np.linspace(0, 100, 10)
            Path(tmp_path / "outputs").mkdir(parents=True, exist_ok=True)

            out[f"slope_{sidx}"][f"dep_{didx}"] = Struct(
                {
                    "waveform_000_deg": Array(wfs),
                    "drift_time_000_deg": Array(drift_time),
                    "drift_time_045_deg": Array(drift_time),
                    "waveform_045_deg": Array(wfs),
                    "r": Array(r),
                    "z": Array(z),
                    "dt": Scalar(1.0),
                }
            )

    info = {}
    info["slope_min"] = Scalar(-1.0)
    info["slope_step"] = Scalar(2.0 / 3)
    info["dep_min"] = Scalar(500)
    info["dep_step"] = Scalar(300 / 3)

    output = {"psl_scan": Struct(out), "info": Struct(info)}

    lh5.write(Struct(output), "V05261B", ideal_psl, wo_mode="of")

    return ideal_psl


@fixture
def make_elecmod(tmp_path):
    out = {"best_fit": {"sigma": 10, "tau": 50}}
    dbetto.utils.write_dict(out, tmp_path / "elecmod.yaml")
    return tmp_path / "elecmod.yaml"


# @mark.needs_remage
def test_drift_time_cli(
    tmp_path,
    monkeypatch,
    legend_stp_path,
    legend_gdml_path,
    make_psl_scan,
    make_elecmod,
):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "extract-drift-time-psl-tuning",
            "--stp-files",
            str(legend_stp_path),
            "--drift-time-file",
            str(tmp_path / "outputs" / "V03422A_drift_times.lh5"),
            "--hpge-detector",
            "V05261B",
            "--psl-file",
            str(make_psl_scan),
            "--elecmod",
            str(make_elecmod),
            "--geom-file",
            str(legend_gdml_path),
        ],
    )

    extract_drift_time_psl_tuning.main()
