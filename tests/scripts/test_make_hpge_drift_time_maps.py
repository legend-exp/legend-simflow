from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest
from lgdo import lh5

testprod = Path(__file__).parent.parent / "dummyprod"
repo_root = Path(__file__).parent.parent.parent


@pytest.mark.needs_julia
@pytest.mark.skipif(shutil.which("julia") is None, reason="julia not installed")
def test_make_hpge_drift_time_maps_l1000(tmp_path):
    """Run the Julia drift time map script and verify the output LH5 structure."""
    dtmap_file = tmp_path / "V05261B-4200V-hpge-drift-time-map.lh5"

    subprocess.run(
        [
            "julia",
            "--project=" + str(repo_root / "workflow/src/LegendSimflow.jl"),
            "--threads",
            "1",
            str(
                repo_root
                / "workflow/src/legendsimflow/scripts/make_hpge_drift_time_maps.jl"
            ),
            "--detector",
            "V05261B",
            "--metadata",
            str(testprod / "inputs"),
            "--opv",
            "4200",
            "--dtmap-settings",
            str(
                testprod
                / "inputs/simprod/config/pars/l1000dsg01/geds/dtmap/settings.yaml"
            ),
            "--output-file",
            str(dtmap_file),
        ],
        check=True,
        cwd=repo_root,
    )

    assert dtmap_file.exists(), "Drift time map LH5 file was not created"

    top_keys = lh5.ls(dtmap_file)
    assert "V05261B" in top_keys, f"Expected group 'V05261B' in LH5, got: {top_keys}"

    inner_keys = lh5.ls(dtmap_file, "V05261B/")
    for expected in ("drift_time_000_deg", "drift_time_045_deg", "r", "z"):
        assert f"V05261B/{expected}" in inner_keys, (
            f"Expected dataset 'V05261B/{expected}' in LH5, found: {inner_keys}"
        )
