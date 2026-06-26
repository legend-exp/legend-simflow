from __future__ import annotations

import sys
from pathlib import Path

import lh5
import numpy as np
import pytest
import yaml
from legendmeta import LegendMetadata
from lgdo import Array, Scalar, Struct
from scipy.stats import norm

from legendsimflow.scripts import (
    build_superpulses_from_data,
    extract_hpge_elec_response_model,
)
from legendsimflow.superpulses import lookup_superpulse_inputs

l200data = Path(__file__).parent.parent / "l200data" / "v3.0.0"


def test_lookup_inputs(test_make_ssc_data):
    meta = LegendMetadata(test_make_ssc_data / "inputs")
    raw_files, evt_files, dsp_config, tab_map, _ = lookup_superpulse_inputs(
        l200data, meta, "l200-p16-r008-ssc", "V03422A", evt_tier_name="pet"
    )

    assert len(raw_files) == 1
    assert len(evt_files) == 1

    assert raw_files[0].exists()
    assert evt_files[0].exists()
    assert dsp_config.exists()

    assert raw_files[0].name == "l200-p16-r008-ssc-20230322T170202Z-tier_raw.lh5"
    assert evt_files[0].name == "l200-p16-r008-ssc-20230322T170202Z-tier_evt.lh5"

    assert tab_map["V03422A"] == 1108804


@pytest.fixture
def test_superpulse_cli(test_make_ssc_data, tmp_path, monkeypatch):
    config_path = tmp_path / "simflow-config.yaml"
    raw_cfg = yaml.safe_load((test_make_ssc_data / "simflow-config.yaml").read_text())
    raw_cfg["paths"]["metadata"] = str(test_make_ssc_data / "inputs")
    raw_cfg["paths"]["l200data"] = str(l200data)
    config_path.write_text(yaml.safe_dump(raw_cfg))

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "build-superpulses-from-data",
            "--runid",
            "l200-p16-r008-ssc",
            "--detector",
            "V03422A",
            "--output-file",
            str(tmp_path / "outputs" / "V03422A_superpulses.lh5"),
            "--plot-file",
            str(tmp_path / "outputs" / "V03422A_superpulses.pdf"),
            "--simflow-config",
            str(config_path),
        ],
    )
    build_superpulses_from_data.main()
    Path(tmp_path / "outputs").mkdir(parents=True, exist_ok=True)
    return str(tmp_path / "outputs" / "V03422A_superpulses.lh5")


@pytest.fixture
def test_make_ideal_psl(tmp_path):
    t = np.arange(5000)
    wfs = []
    # Gaussian PDF
    for mu in np.linspace(0, 1990, 200):
        wfs.append(np.cumsum(norm.pdf(t, loc=mu, scale=100)))
        wfs[-1] /= wfs[-1][-1]  # normalize to 1 at the end of the waveform

    wfs = np.array(wfs).reshape(20, 10, 5000)
    Path(tmp_path / "outputs").mkdir(parents=True, exist_ok=True)

    ideal_psl = str(tmp_path / "outputs" / "l200-p16-r008-ssc-ideal_psl.lh5")

    output = Struct({"waveform_000_deg": Array(wfs), "dt": Scalar(1.0)})
    lh5.write(output, "V03422A", ideal_psl, wo_mode="of")

    return ideal_psl


def test_extract_electronics_model_cli_with_data(
    test_make_ssc_data, tmp_path, monkeypatch, test_make_ideal_psl, test_superpulse_cli
):
    # test first with defaults (legend)

    config_path = tmp_path / "simflow-config.yaml"
    raw_cfg = yaml.safe_load((test_make_ssc_data / "simflow-config.yaml").read_text())
    raw_cfg["paths"]["metadata"] = str(test_make_ssc_data / "inputs")
    config_path.write_text(yaml.safe_dump(raw_cfg))

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "extract-hpge-electronics-model",
            "--runid",
            "l200-p16-r008-ssc",
            "--hpge-detector",
            "V03422A",
            "--pars-file",
            str(tmp_path / "outputs" / "V03422A_electronics_pars.yaml"),
            "--simflow-config",
            str(config_path),
        ],
    )

    extract_hpge_elec_response_model.main()

    pars_file = tmp_path / "outputs" / "V03422A_electronics_pars.yaml"
    assert pars_file.exists()
    pars = yaml.safe_load(pars_file.read_text())

    assert "sigma" in pars
    assert "tau" in pars

    assert isinstance(pars["sigma"], (int, float))
    assert isinstance(pars["tau"], (int, float))

    # now test with non-default settings (l200)

    config_path = tmp_path / "simflow-config.yaml"
    raw_cfg = yaml.safe_load((test_make_ssc_data / "simflow-config.yaml").read_text())

    raw_cfg["paths"]["metadata"] = str(test_make_ssc_data / "inputs")
    raw_cfg["paths"]["l200data"] = str(l200data)
    raw_cfg["experiment"] = "l200cfg01"

    config_path.write_text(yaml.safe_dump(raw_cfg))

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "extract-hpge-electronics-model",
            "--runid",
            "l200-p16-r008-ssc",
            "--hpge-detector",
            "V03422A",
            "--ideal-lib",
            test_make_ideal_psl,
            "--superpulses",
            test_superpulse_cli,
            "--pars-file",
            str(tmp_path / "outputs" / "V03422A_electronics_pars.yaml"),
            "--plot-file",
            str(tmp_path / "outputs" / "V03422A_electronics_tuning.pdf"),
            "--simflow-config",
            str(config_path),
        ],
    )

    extract_hpge_elec_response_model.main()

    pars_file = tmp_path / "outputs" / "V03422A_electronics_pars.yaml"
    assert pars_file.exists()
    pars = yaml.safe_load(pars_file.read_text())

    assert "sigma" in pars
    assert "tau" in pars

    assert isinstance(pars["sigma"], (int, float))
    assert isinstance(pars["tau"], (int, float))

    plot_file = tmp_path / "outputs" / "V03422A_electronics_tuning.pdf"
    print(plot_file)
    assert plot_file.exists()
