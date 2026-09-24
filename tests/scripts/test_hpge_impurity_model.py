from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from lgdo import Array
import lh5
import numpy as np
import pytest

from legendsimflow.psl import validate_ssd_scan_grid
from legendsimflow.superpulses import lookup_superpulse_inputs
from legendmeta import LegendMetadata
testprod = Path(__file__).parent.parent / "dummyprod"
repo_root = Path(__file__).parent.parent.parent

l200data = Path(__file__).parent.parent / "l200data" / "v3.0.0"
DETECTOR = "V03422A"


def _drift_time(peak1,peak2, frac = 0.5, size = 1000):
    """Get a drift time distribution from two peaks with a given fraction of events in each peak."""
    n1 = int(size * frac)
    n2 = size - n1
    peak1 = np.random.normal(peak1, 5, n1)
    peak2 = np.random.normal(peak2, 5, n2)
    return np.concatenate([peak1, peak2])

@pytest.fixture
def make_sim_drift_time(tmp_path):
    
    dt = Array(_drift_time(peak1=1000, peak2=2000, frac=0.5, size=10000))
    sim_file = tmp_path / "sim.lh5"
    lh5.write(dt,"drift_time", sim_file)

    return sim_file

def test_hpge_impurity_cli_with_data(test_make_ssc_data, make_sim_drift_time, tmp_path):
    meta = LegendMetadata(test_make_ssc_data / "inputs", lazy = True)
    _, evt_files, _, _, _ = lookup_superpulse_inputs(
        l200data, meta, "l200-p16-r008-ssc", DETECTOR, evt_tier_name="pet"
    )

    lh5.show(evt_files[0])

    print()

