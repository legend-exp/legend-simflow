from __future__ import annotations

import sys
from pathlib import Path

import dbetto
import lh5
import numpy as np
import pytest
from legendmeta import LegendMetadata
from lgdo import Array, Scalar, Struct

from legendsimflow import utils
from legendsimflow.impurity_tuning import get_run_mapping
from legendsimflow.metadata import get_simconfig
from legendsimflow.scripts.extract_hpge_impurity_model import main
from legendsimflow.superpulses import lookup_superpulse_inputs

testprod = Path(__file__).parent.parent / "dummyprod"
repo_root = Path(__file__).parent.parent.parent

l200data = Path(__file__).parent.parent / "l200data" / "v3.0.0"
DETECTOR = "V03422A"
rng = np.random.default_rng(seed=42)


def _drift_time(peak1, peak2, frac=0.5, size=1000):
    """Get a drift time distribution from two peaks with a given fraction of events in each peak."""
    n1 = int(size * frac)
    n2 = size - n1
    peak1 = rng.normal(peak1, 100, n1)
    peak2 = rng.normal(peak2, 100, n2)
    return np.concatenate([peak1, peak2])


def test_get_run_mapping(test_make_ssc_data):
    config = utils.init_simflow_context(
        test_make_ssc_data / "simflow-config-l200-ssc.yaml", workflow=None
    ).config
    runs = ["r008"]
    mapping = get_run_mapping(get_simconfig(config, "hit", simid=None), runs)

    assert mapping == {"l200-p16-r008-ssc": "source_pos_1"}


@pytest.fixture
def make_sim_drift_time(tmp_path):
    sim_file = tmp_path / "source_pos_1.lh5"

    for slope in range(1, 10):
        for depv in range(1, 10):
            dt = Array(
                _drift_time(
                    peak1=1000 + 10 * (slope - 5),
                    peak2=2000 + 10 * (depv - 5),
                    frac=0.5,
                    size=10000,
                )
            )

            wo_mode = "of" if slope == 1 and depv == 1 else "append_column"
            lh5.write(
                dt,
                f"V03422A/psl_scan/slope_{slope}/dep_{depv}/drift_time",
                sim_file,
                wo_mode=wo_mode,
            )

    lh5.write(
        Array(rng.uniform(500, 3000, size=10000)),
        "V03422A/energy",
        sim_file,
        wo_mode="append_column",
    )

    grid_info = Struct(
        {
            "slope_min": Scalar(-1),
            "slope_step": Scalar(0.1),
            "dep_min": Scalar(1000),
            "dep_step": Scalar(100),
        }
    )
    lh5.write(grid_info, "V03422A/grid_info", sim_file, wo_mode="append_column")
    return sim_file


def test_hpge_impurity_cli_with_data(
    test_make_ssc_data, make_sim_drift_time, tmp_path, monkeypatch
):
    meta = LegendMetadata(test_make_ssc_data / "inputs", lazy=True)
    _, evt_files, _, _, _ = lookup_superpulse_inputs(
        l200data, meta, "l200-p16-r008-ssc", DETECTOR, evt_tier_name="pet"
    )

    run_norms = {"l200-p16-r008-ssc": 1.0}
    Path(tmp_path / "outputs").mkdir(parents=True)
    dbetto.utils.write_dict(run_norms, tmp_path / "outputs" / "run_norms.yaml")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "extract-hpge-impurity-model",
            "--drift-time",
            str(make_sim_drift_time),
            "--pars-file",
            str(tmp_path / "outputs" / f"{DETECTOR}_electronics_pars.yaml"),
            "--simflow-config",
            str(test_make_ssc_data / "simflow-config-l200-ssc.yaml"),
            "--plot-file",
            str(tmp_path / "outputs" / f"{DETECTOR}_drift_time_obs.pdf"),
            "--data-path",
            str(evt_files[0].parent),
            "--pars-file",
            str(tmp_path / "outputs" / f"{DETECTOR}_electronics_pars.yaml"),
            "--simflow-config",
            str(test_make_ssc_data / "simflow-config-l200-ssc.yaml"),
            "--runids",
            "l200-p16-r008-ssc",
        ],
    )

    main()

    assert (tmp_path / "outputs" / f"{DETECTOR}_electronics_pars.yaml").exists()
    assert (tmp_path / "outputs" / f"{DETECTOR}_drift_time_obs.pdf").exists()
