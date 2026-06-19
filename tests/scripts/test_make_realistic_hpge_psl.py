from __future__ import annotations

import sys
from pathlib import Path

import lh5
import numpy as np
import yaml
from lgdo import Array, Scalar, Struct

from legendsimflow.scripts import make_hpge_realistic_pulse_shape_lib

dummyprod = Path(__file__).parent.parent / "dummyprod"


def test_make_realistic_hpge_psl_cli(tmp_path, monkeypatch):
    detector = "V05261B"

    config_path = tmp_path / "simflow-config.yaml"
    raw_cfg = yaml.safe_load((dummyprod / "simflow-config.yaml").read_text())
    raw_cfg["paths"]["metadata"] = str(dummyprod / "inputs")
    config_path.write_text(yaml.safe_dump(raw_cfg))

    t = np.arange(5000)
    wfs = []
    for t0 in np.linspace(400, 2000, 6):
        wf = np.clip((t - t0) / 250.0, 0.0, 1.0)
        wfs.append(wf)
    wfs = np.asarray(wfs).reshape(2, 3, 5000)

    ideal_psl_file = tmp_path / "ideal_psl.lh5"
    ideal_struct = Struct(
        {
            "r": Array(np.array([0.01, 0.02]), attrs={"units": "m"}),
            "z": Array(np.array([-0.01, 0.00, 0.01]), attrs={"units": "m"}),
            "waveform_000_deg": Array(wfs),
            "dt": Scalar(1.0, attrs={"units": "ns"}),
        }
    )
    lh5.write(ideal_struct, detector, ideal_psl_file, wo_mode="of")

    electronics_model_file = tmp_path / "electronics_model.yaml"
    electronics_model_file.write_text(
        yaml.safe_dump({detector: {"sigma": 60.0, "tau": 100.0}})
    )

    output_file = tmp_path / "realistic_psl.lh5"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "make-hpge-realistic-pulse-shape-lib",
            "--detector",
            detector,
            "--electronics-model-file",
            str(electronics_model_file),
            "--plot-file",
            str(tmp_path / "validation_plots.pdf"),
            "--input-file",
            str(ideal_psl_file),
            "--output-file",
            str(output_file),
            "--simflow-config",
            str(config_path),
        ],
    )
    print(str(tmp_path / "validation_plots.pdf"))

    make_hpge_realistic_pulse_shape_lib.main()

    assert output_file.exists(), "Realistic PSL output file was not created"

    top_keys = lh5.ls(output_file)
    assert detector in top_keys

    out_struct = lh5.read(detector, output_file)
    expected_keys = {"r", "z", "dt", "t0", "waveform_000_deg", "drift_time_000_deg"}
    assert expected_keys.issubset(set(out_struct.keys()))

    assert out_struct["waveform_000_deg"].view_as("np").shape == (2, 3, 4001)
    assert out_struct["drift_time_000_deg"].view_as("np").shape == (2, 3)
