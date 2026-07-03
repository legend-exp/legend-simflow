from __future__ import annotations

import sys
from pathlib import Path

import lh5
import numpy as np
import pytest
import yaml

from legendsimflow.scripts.tier import evt

dummyprod = Path(__file__).parent.parent / "dummyprod"


@pytest.mark.needs_remage
def test_evt_script_cli(
    tmp_path,
    monkeypatch,
    legend_stp_path,
    legend_opt_path,
    legend_hit_path,
    legend_simstat_part_path,
    legend_detector_usabilities_path,
):
    raw = yaml.safe_load((dummyprod / "simflow-config-l1000.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    config_path = tmp_path / "simflow-config-l1000.yaml"
    config_path.write_text(yaml.safe_dump(raw))

    evt_file = tmp_path / "evt.lh5"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "evt",
            "--stp-file",
            str(legend_stp_path),
            "--opt-file",
            str(legend_opt_path),
            "--hit-file",
            str(legend_hit_path),
            "--simstat-part-file",
            str(legend_simstat_part_path),
            "--detector-usabilities-file",
            str(legend_detector_usabilities_path),
            "--jobid",
            "0000",
            "--evt-file",
            str(evt_file),
            "--simflow-config",
            str(config_path),
        ],
    )
    evt.main()

    assert evt_file.exists(), "evt output file was not created"

    root_keys = lh5.ls(evt_file)
    assert "evt" in root_keys, f"'evt' group missing from {evt_file}; got {root_keys}"

    evt_groups = lh5.ls(evt_file, "evt/")
    for group in ("evt/trigger", "evt/geds", "evt/spms", "evt/coincident"):
        assert group in evt_groups, f"'{group}' missing from evt/; got {evt_groups}"

    trigger_fields = {
        f.removeprefix("evt/trigger/") for f in lh5.ls(evt_file, "evt/trigger/")
    }
    for field in ("run", "period", "evtid", "timestamp"):
        assert field in trigger_fields, (
            f"'trigger/{field}' missing; got {trigger_fields}"
        )

    geds_fields = {f.removeprefix("evt/geds/") for f in lh5.ls(evt_file, "evt/geds/")}
    for field in (
        "energy",
        "energy_sum",
        "is_good_channel",
        "multiplicity",
        "rawid",
        "hit_idx",
        "psd",
    ):
        assert field in geds_fields, f"'geds/{field}' missing; got {geds_fields}"

    psd_fields = {
        f.removeprefix("evt/geds/psd/") for f in lh5.ls(evt_file, "evt/geds/psd/")
    }
    for field in ("is_good", "is_valid_sim", "aoe", "has_aoe", "is_single_site"):
        assert field in psd_fields, f"'geds/psd/{field}' missing; got {psd_fields}"

    spms_fields = {f.removeprefix("evt/spms/") for f in lh5.ls(evt_file, "evt/spms/")}
    for field in (
        "energy",
        "energy_sum",
        "is_saturated",
        "rawid",
        "hit_idx",
        "time",
        "multiplicity",
    ):
        assert field in spms_fields, f"'spms/{field}' missing; got {spms_fields}"

    coincident_fields = {
        f.removeprefix("evt/coincident/") for f in lh5.ls(evt_file, "evt/coincident/")
    }
    for field in ("geds", "spms"):
        assert field in coincident_fields, (
            f"'coincident/{field}' missing; got {coincident_fields}"
        )

    # RC is disabled → rc_energy and rc_time must not be present
    assert "rc_energy" not in spms_fields, (
        "spms/rc_energy present but --add-random-coincidences was not passed"
    )
    assert "rc_time" not in spms_fields, (
        "spms/rc_time present but --add-random-coincidences was not passed"
    )

    def _read(path: str) -> np.ndarray:
        return lh5.read_as(path, str(evt_file), library="np")

    # geds/energy_sum: float32, all finite and non-negative
    geds_esum = _read("evt/geds/energy_sum")
    assert geds_esum.dtype == np.float32, (
        f"geds/energy_sum dtype should be float32, got {geds_esum.dtype}"
    )
    assert np.all(np.isfinite(geds_esum)), "geds/energy_sum contains non-finite values"
    assert np.all(geds_esum >= 0), "geds/energy_sum contains negative values"

    # spms/energy_sum: float32, all finite and non-negative
    spms_esum = _read("evt/spms/energy_sum")
    assert spms_esum.dtype == np.float32, (
        f"spms/energy_sum dtype should be float32, got {spms_esum.dtype}"
    )
    assert np.all(np.isfinite(spms_esum)), "spms/energy_sum contains non-finite values"
    assert np.all(spms_esum >= 0), "spms/energy_sum contains negative values"

    # geds/multiplicity: non-negative integers
    geds_mult = _read("evt/geds/multiplicity")
    assert np.all(geds_mult >= 0), "geds/multiplicity contains negative values"

    # spms/multiplicity: non-negative integers
    spms_mult = _read("evt/spms/multiplicity")
    assert np.all(spms_mult >= 0), "spms/multiplicity contains negative values"

    # trigger/period: all p03 → code 3
    period_vals = _read("evt/trigger/period")
    assert np.all(period_vals == 3), (
        f"expected period=3 (p03), got unique={np.unique(period_vals)}"
    )

    # trigger/run: values in {0, 1} for r000 and r001
    run_vals = _read("evt/trigger/run")
    assert np.all(np.isin(run_vals, [0, 1])), (
        f"unexpected run values: {np.unique(run_vals)}"
    )


@pytest.mark.needs_remage
def test_evt_script_cli_skip_opt(
    tmp_path,
    monkeypatch,
    legend_stp_path,
    legend_hit_path,
    legend_simstat_part_path,
    legend_detector_usabilities_path,
):
    """--skip-opt: no opt file passed; geds table present, spms table absent."""
    raw = yaml.safe_load((dummyprod / "simflow-config-l1000.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    config_path = tmp_path / "simflow-config-l1000.yaml"
    config_path.write_text(yaml.safe_dump(raw))

    evt_file = tmp_path / "evt.lh5"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "evt",
            "--stp-file",
            str(legend_stp_path),
            "--hit-file",
            str(legend_hit_path),
            "--simstat-part-file",
            str(legend_simstat_part_path),
            "--detector-usabilities-file",
            str(legend_detector_usabilities_path),
            "--jobid",
            "0000",
            "--evt-file",
            str(evt_file),
            "--simflow-config",
            str(config_path),
            "--skip-opt",
        ],
    )
    evt.main()

    assert evt_file.exists(), "evt output file was not created"

    evt_groups = lh5.ls(evt_file, "evt/")

    # these groups must be present
    for group in ("evt/trigger", "evt/geds", "evt/coincident"):
        assert group in evt_groups, f"'{group}' missing from evt/; got {evt_groups}"

    # spms table must be absent (opt tier was skipped)
    assert "evt/spms" not in evt_groups, (
        f"'evt/spms' should be absent when --skip-opt is set, got {evt_groups}"
    )

    coincident_fields = {
        f.removeprefix("evt/coincident/") for f in lh5.ls(evt_file, "evt/coincident/")
    }

    # geds coincidence flag must be present
    assert "geds" in coincident_fields, (
        f"'coincident/geds' missing; got {coincident_fields}"
    )
    # spms coincidence flag must be absent (no opt data)
    assert "spms" not in coincident_fields, (
        f"'coincident/spms' should be absent when --skip-opt is set; got {coincident_fields}"
    )


@pytest.mark.needs_remage
def test_evt_script_cli_skip_hit(
    tmp_path,
    monkeypatch,
    legend_stp_path,
    legend_opt_path,
    legend_simstat_part_path,
    legend_detector_usabilities_path,
):
    """--skip-hit: no hit file passed; spms table present, geds table absent.

    The trigger fields must be sourced from the opt tier, so trigger/period
    should still be 3 (p03 runs).
    """
    raw = yaml.safe_load((dummyprod / "simflow-config-l1000.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    config_path = tmp_path / "simflow-config-l1000.yaml"
    config_path.write_text(yaml.safe_dump(raw))

    evt_file = tmp_path / "evt.lh5"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "evt",
            "--stp-file",
            str(legend_stp_path),
            "--opt-file",
            str(legend_opt_path),
            "--simstat-part-file",
            str(legend_simstat_part_path),
            "--detector-usabilities-file",
            str(legend_detector_usabilities_path),
            "--jobid",
            "0000",
            "--evt-file",
            str(evt_file),
            "--simflow-config",
            str(config_path),
            "--skip-hit",
        ],
    )
    evt.main()

    assert evt_file.exists(), "evt output file was not created"

    evt_groups = lh5.ls(evt_file, "evt/")

    # these groups must be present
    for group in ("evt/trigger", "evt/spms", "evt/coincident"):
        assert group in evt_groups, f"'{group}' missing from evt/; got {evt_groups}"

    # geds table must be absent (hit tier was skipped)
    assert "evt/geds" not in evt_groups, (
        f"'evt/geds' should be absent when --skip-hit is set, got {evt_groups}"
    )

    coincident_fields = {
        f.removeprefix("evt/coincident/") for f in lh5.ls(evt_file, "evt/coincident/")
    }

    # spms coincidence flag must be present
    assert "spms" in coincident_fields, (
        f"'coincident/spms' missing; got {coincident_fields}"
    )
    # geds coincidence flag must be absent (no hit data)
    assert "geds" not in coincident_fields, (
        f"'coincident/geds' should be absent when --skip-hit is set; got {coincident_fields}"
    )

    # with --skip-hit the trigger is sourced from the opt tier; the opt fixture
    # uses l200-p03 runs so period must always be 3
    period_vals = lh5.read_as("evt/trigger/period", str(evt_file), library="np")
    assert np.all(period_vals == 3), (
        f"expected period=3 (p03, from opt tier), got unique={np.unique(period_vals)}"
    )


def test_evt_script_cli_skip_both_raises(monkeypatch):
    """Passing both --skip-opt and --skip-hit must cause an immediate SystemExit."""
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "evt",
            "--stp-file",
            "dummy_stp.lh5",
            "--simstat-part-file",
            "dummy_simstat.yaml",
            "--detector-usabilities-file",
            "dummy_usabilities.yaml",
            "--jobid",
            "0000",
            "--evt-file",
            "dummy_evt.lh5",
            "--simflow-config",
            "dummy_config.yaml",
            "--skip-opt",
            "--skip-hit",
        ],
    )
    with pytest.raises(SystemExit):
        evt.main()
