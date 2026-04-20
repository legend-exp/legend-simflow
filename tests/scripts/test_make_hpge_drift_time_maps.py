from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest
from lgdo import lh5

# Repo root is one level above the tests/ directory.
repo_root = Path(__file__).parent.parent.parent

dummyprod_inputs = Path(__file__).parent.parent / "dummyprod" / "inputs"

julia_script = (
    repo_root
    / "workflow"
    / "src"
    / "legendsimflow"
    / "scripts"
    / "make_hpge_drift_time_maps.jl"
)
julia_project = repo_root / "workflow" / "src" / "LegendSimflow.jl"


@pytest.mark.needs_julia
@pytest.mark.skipif(shutil.which("julia") is None, reason="julia not installed")
def test_make_hpge_drift_time_maps_l1000(legend_dtmap_path):
    assert legend_dtmap_path.exists(), "Drift time map LH5 file was not created"

    top_keys = lh5.ls(legend_dtmap_path)
    assert "V05261B" in top_keys, f"Expected group 'V05261B' in LH5, got: {top_keys}"

    inner_keys = lh5.ls(legend_dtmap_path, "V05261B/")
    for expected in ("drift_time_000_deg", "drift_time_045_deg", "r", "z"):
        assert f"V05261B/{expected}" in inner_keys, (
            f"Expected dataset 'V05261B/{expected}' in LH5, found: {inner_keys}"
        )


@pytest.mark.needs_julia
@pytest.mark.skipif(shutil.which("julia") is None, reason="julia not installed")
def test_make_hpge_drift_time_maps_produces_valid_lh5(tmp_path):
    output_file = tmp_path / "dtmap_V99000A.lh5"

    cmd = [
        "julia",
        "--project=" + str(julia_project),
        "--threads",
        "1",
        str(julia_script),
        "--detector",
        "V99000A",
        "--metadata",
        str(dummyprod_inputs),
        "--opv",
        "4200",
        "--dtmap-settings",
        str(dummyprod_inputs / "simprod/config/pars/legend/geds/dtmap/settings.yaml"),
        "--output-file",
        str(output_file),
    ]

    result = subprocess.run(
        cmd,
        check=False,
        cwd=repo_root,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        pytest.fail(
            f"Julia script exited with code {result.returncode}.\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )

    assert output_file.exists(), "Output LH5 file was not created"

    top_keys = lh5.ls(output_file)
    assert "V99000A" in top_keys, f"Expected group 'V99000A' in LH5, got: {top_keys}"

    inner_keys = lh5.ls(output_file, "V99000A/")
    for expected in ("drift_time_000_deg", "drift_time_045_deg", "r", "z"):
        assert f"V99000A/{expected}" in inner_keys, (
            f"Expected dataset 'V99000A/{expected}' in LH5, found: {inner_keys}"
        )
