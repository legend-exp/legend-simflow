from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import lh5
import pytest

testprod = Path(__file__).parent.parent / "dummyprod"
repo_root = Path(__file__).parent.parent.parent


@pytest.mark.needs_julia
@pytest.mark.skipif(shutil.which("julia") is None, reason="julia not installed")
def test_make_hpge_pulse_shape_lib_l1000(tmp_path):
    """Run the Julia drift time map script and verify the output LH5 structure."""
    psl_file = tmp_path / "V05261B-4200V-hpge-pulse-shape-lib.lh5"

    subprocess.run(
        [
            "julia",
            "--project=" + str(repo_root / "workflow/src/LegendSimflow.jl"),
            "--threads",
            "1",
            str(
                repo_root
                / "workflow/src/legendsimflow/scripts/make_hpge_ideal_pulse_shape_lib.jl"
            ),
            "--detector",
            "V05261B",
            "--metadata",
            str(testprod / "inputs"),
            "--opv",
            "4200",
            "--ssd-settings",
            str(
                testprod
                / "inputs/simprod/config/pars/l1000dsg01/geds/ssd/settings.yaml"
            ),
            "--output-file",
            str(psl_file),
        ],
        check=True,
        cwd=repo_root,
    )

    assert psl_file.exists(), "Pulse shape library LH5 file was not created"

    top_keys = lh5.ls(psl_file)
    assert "V05261B" in top_keys, f"Expected group 'V05261B' in LH5, got: {top_keys}"
