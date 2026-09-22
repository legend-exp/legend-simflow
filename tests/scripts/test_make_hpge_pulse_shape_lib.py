from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import lh5
import numpy as np
import pytest

from legendsimflow.psl import validate_ssd_scan_grid

testprod = Path(__file__).parent.parent / "dummyprod"
repo_root = Path(__file__).parent.parent.parent


@pytest.mark.needs_julia
@pytest.mark.skipif(shutil.which("julia") is None, reason="julia not installed")
def test_make_hpge_pulse_shape_lib_l200(tmp_path):
    psl_file = tmp_path / "V00001A-3500V-hpge-pulse-shape-lib.lh5"

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
            "V00001A",
            "--metadata",
            str(testprod / "legend-metadata"),
            "--opv",
            "3500",
            "--ssd-settings",
            str(
                testprod
                / "legend-metadata/simprod/config/pars/l200cfg01/geds/ssd/settings.yaml"
            ),
            "--output-file",
            str(psl_file),
        ],
        check=True,
        cwd=repo_root,
    )

    assert psl_file.exists(), "Pulse shape library LH5 file was not created"

    top_keys = lh5.ls(psl_file)
    assert "V00001A" in top_keys, f"Expected group 'V00001A' in LH5, got: {top_keys}"


@pytest.mark.needs_julia
@pytest.mark.skipif(shutil.which("julia") is None, reason="julia not installed")
def test_make_hpge_pulse_shape_lib_scan(tmp_path):
    psl_file = tmp_path / "V05261B-4200V-hpge-pulse-shape-scan-lib.lh5"

    subprocess.run(
        [
            "julia",
            "--project=" + str(repo_root / "workflow/src/LegendSimflow.jl"),
            "--threads",
            "1",
            str(
                repo_root
                / "workflow/src/legendsimflow/scripts/make_hpge_ideal_pulse_shape_lib_scan.jl"
            ),
            "--detector",
            "V05261B",
            "--metadata",
            str(testprod / "inputs"),
            "--opv",
            "4200",
            "--ssd-settings",
            str(testprod / "inputs/simprod/config/pars/legend/geds/ssd/settings.yaml"),
            "--scan-settings",
            str(
                testprod
                / "inputs/simprod/config/pars/legend/geds/ssd/scan_settings.yaml"
            ),
            "--output-file",
            str(psl_file),
        ],
        check=True,
        cwd=repo_root,
    )

    assert validate_ssd_scan_grid(str(psl_file), "V05261B")
    assert psl_file.exists(), "Pulse shape library LH5 file was not created"

    top_keys = lh5.ls(psl_file)

    assert "V05261B" in top_keys, f"Expected group 'V05261B' in LH5, got: {top_keys}"

    def names(group):
        return {key.rsplit("/", 1)[-1] for key in lh5.ls(psl_file, group)}

    # the scan grid is pinned by scan_settings.yaml, slope "-1:1:0" and
    # depv_shift "-500:450:-50", so both dimensions are two points wide
    assert names("V05261B/") == {"psl_scan", "grid_info"}
    assert names("V05261B/psl_scan/") == {"slope_1", "slope_2"}
    assert names("V05261B/psl_scan/slope_1/") == {"dep_1", "dep_2"}
    assert names("V05261B/psl_scan/slope_1/dep_1/") == {
        "r",
        "z",
        "dt",
        "waveform_000_deg",
        "waveform_045_deg",
        "impurity_scale",
    }

    def waveform(slope, dep):
        path = f"V05261B/psl_scan/slope_{slope}/dep_{dep}/waveform_000_deg"
        return np.asarray(lh5.read(path, psl_file))

    # the scan must move the simulation along both axes, otherwise it is
    # writing the same detector over and over. equal_nan is required: a quarter
    # of every map is the NaN padding outside the detector, and without it
    # array_equal is False even for a map compared against itself
    assert not np.array_equal(waveform(1, 1), waveform(2, 1), equal_nan=True), (
        "the impurity slope does not change the waveforms"
    )
    assert not np.array_equal(waveform(1, 1), waveform(1, 2), equal_nan=True), (
        "the depletion voltage does not change the waveforms"
    )

    # info holds the first point of each range, the depletion voltage as an
    # absolute value rather than a shift
    assert lh5.read("V05261B/grid_info/dep_min", psl_file).value == 4200 - 500
    assert lh5.read("V05261B/grid_info/slope_min", psl_file).value == -1
