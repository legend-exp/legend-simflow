from __future__ import annotations

import sys
from pathlib import Path

import lh5
import numpy as np
import pytest
import yaml

from legendsimflow.scripts.tier import opt

dummyprod = Path(__file__).parent.parent / "dummyprod"

# channels with non-trivial usabilities added to the dummyprod status file
_SIPM_ON = "S002"  # usability: on  → code 0
_SIPM_AC = "S003"  # usability: ac  → code 1
_SIPM_OFF = "S007"  # usability: off → code 2


@pytest.mark.needs_remage
def test_opt_script_cli(
    tmp_path,
    monkeypatch,
    legend_testdata,
    legend_stp_path,
    legend_gdml_path,
    legend_simstat_part_path,
    legend_detector_usabilities_path,
):
    optmap_path = Path(legend_testdata.get_path("remage/l200cfg01-optmap-dummy.lh5"))

    raw = yaml.safe_load((dummyprod / "simflow-config-l1000.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    config_path = tmp_path / "simflow-config-l1000.yaml"
    config_path.write_text(yaml.safe_dump(raw))

    common_args = [
        "--stp-file", str(legend_stp_path),
        "--optmap-lar", str(optmap_path),
        "--geom-file", str(legend_gdml_path),
        "--detector-usabilities-file", str(legend_detector_usabilities_path),
        "--jobid", "0000",
        "--scintillator-volume-name", "liquid_argon",
        "--simflow-config", str(config_path),
    ]  # fmt: skip

    # --- run 1: sum optical map (all SiPMs together) ---
    opt_file = tmp_path / "opt.lh5"
    monkeypatch.setattr(sys, "argv", [
        "opt", *common_args,
        "--simstat-part-file", str(legend_simstat_part_path),
        "--opt-file", str(opt_file),
    ])  # fmt: skip
    opt.main()

    assert opt_file.exists(), "opt output file was not created"

    root_keys = lh5.ls(opt_file)
    assert "hit" in root_keys, f"'hit' group missing; got {root_keys}"
    assert "tcm" in root_keys, f"'tcm' group missing; got {root_keys}"
    assert "hit/spms" in lh5.ls(opt_file, "hit/"), (
        f"'hit/spms' missing; got {lh5.ls(opt_file, 'hit/')}"
    )

    expected_fields = {
        "energy", "evtid", "is_saturated", "period", "run", "t0", "time", "usability",
    }  # fmt: skip
    spms_fields = {f.removeprefix("hit/spms/") for f in lh5.ls(opt_file, "hit/spms/")}
    assert not (expected_fields - spms_fields), (
        f"fields {expected_fields - spms_fields} missing from hit/spms"
    )

    def _field(table: str, name: str, f: Path = opt_file) -> np.ndarray:
        return lh5.read_as(f"hit/{table}/{name}", f, library="np")

    assert np.all(_field("spms", "period") == 3)
    assert np.all(np.isin(_field("spms", "run"), [0, 1]))
    assert np.all(np.isin(_field("spms", "usability"), [0, 1, 2]))

    energy = _field("spms", "energy")
    assert np.any(np.isfinite(energy)), "energy contains no finite values"
    assert np.all(energy[np.isfinite(energy)] >= 0), "energy contains negative values"

    is_sat = _field("spms", "is_saturated")
    assert is_sat.dtype == np.bool_, f"unexpected is_saturated dtype {is_sat.dtype}"

    # --- run 2: per-SiPM mode with tiny partition → check usability codes ---
    part_file = tmp_path / "partitions_small.yaml"
    part_file.write_text(yaml.safe_dump({"job_0000": {"l200-p03-r000-phy": [0, 99]}}))

    opt_per_sipm = tmp_path / "opt_per_sipm.lh5"
    monkeypatch.setattr(sys, "argv", [
        "opt", *common_args,
        "--simstat-part-file", str(part_file),
        "--opt-file", str(opt_per_sipm),
        "--optmap-per-sipm",
    ])  # fmt: skip
    opt.main()

    hit_tables = {t.removeprefix("hit/") for t in lh5.ls(opt_per_sipm, "hit/")}
    for sipm in (_SIPM_ON, _SIPM_AC, _SIPM_OFF):
        assert sipm in hit_tables, f"expected hit/{sipm} not found; got {hit_tables}"

    assert np.all(_field(_SIPM_ON, "usability", opt_per_sipm) == 0), (
        f"{_SIPM_ON} usability should be 0 (on)"
    )
    assert np.all(_field(_SIPM_AC, "usability", opt_per_sipm) == 1), (
        f"{_SIPM_AC} usability should be 1 (ac)"
    )
    assert np.all(_field(_SIPM_OFF, "usability", opt_per_sipm) == 2), (
        f"{_SIPM_OFF} usability should be 2 (off)"
    )
