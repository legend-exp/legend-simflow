from __future__ import annotations

import sys
from pathlib import Path

import yaml
from legendmeta import LegendMetadata

from legendsimflow.scripts.build_superpulses_from_data import main
from legendsimflow.superpulses import lookup_superpulse_inputs

l200data = Path(__file__).parent.parent / "l200data" / "v3.0.0"


def test_lookup_inputs(test_make_ssc_data):
    meta = LegendMetadata(test_make_ssc_data / "inputs")
    raw_files, evt_files, dsp_config, tab_map = lookup_superpulse_inputs(
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
    main()


def test_extract_electronics_model_cli():
    pass
