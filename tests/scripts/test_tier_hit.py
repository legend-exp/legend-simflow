from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml
from lgdo import lh5

from legendsimflow.scripts.tier import hit

dummyprod = Path(__file__).parent.parent / "dummyprod"

_STP_FILE = (
    dummyprod
    / "generated-l1000/tier/stp/ultem_insulators_Pb214_to_Po214"
    / "l1000dsg01-ultem_insulators_Pb214_to_Po214-job_0000-tier_stp.lh5"
)


@pytest.mark.skipif(
    not _STP_FILE.exists(),
    reason="generated-l1000 outputs not available (run test_l1000_workflow first)",
)
def test_hit_script_cli(tmp_path, monkeypatch):
    # Build a config with all $_ paths resolved to absolute dummyprod paths so
    # that patterns.output_eresmod_filename() etc. point to the pre-generated
    # parameter files inside generated-l1000/.
    generated = dummyprod / "generated-l1000"

    config_path = tmp_path / "simflow-config-l1000.yaml"
    raw = yaml.safe_load((dummyprod / "simflow-config-l1000.yaml").read_text())

    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    raw["paths"]["generated"] = str(generated)
    raw["paths"]["benchmarks"] = str(generated / "benchmarks")
    raw["paths"]["log"] = str(generated / "log")
    raw["paths"]["config"] = str(dummyprod / "inputs/simprod/config")
    raw["paths"]["pars"] = str(generated / "pars")
    raw["paths"]["macros"] = str(generated / "macros")
    raw["paths"]["tier"] = {
        "vtx": str(generated / "tier/vtx"),
        "stp": str(generated / "tier/stp"),
        "opt": str(generated / "tier/opt"),
        "hit": str(generated / "tier/hit"),
        "evt": str(generated / "tier/evt"),
        "cvt": str(generated / "tier/cvt"),
        "pdf": str(generated / "tier/pdf"),
    }
    raw["paths"]["pdf_releases"] = str(generated / "releases")
    raw["paths"]["optical_maps"] = {
        "fiber": str(dummyprod / "inputs/simprod/l200a-optical-map-fiber.lh5"),
        "lar": str(dummyprod / "inputs/simprod/l200a-optical-map-lar.lh5"),
        "pen": str(dummyprod / "inputs/simprod/l200a-optical-map-pen.lh5"),
    }
    config_path.write_text(yaml.safe_dump(raw))

    hit_file = tmp_path / "hit.lh5"

    argv = [
        "hit",
        "--stp-file",
        str(_STP_FILE),
        "--jobid",
        "0000",
        "--hit-file",
        str(hit_file),
        "--geom-file",
        str(
            generated
            / "pars/geom/l1000dsg01-ultem_insulators_Pb214_to_Po214-tier_stp-geom.gdml"
        ),
        "--dtmap-files",
        "--currmod-files",
        "--simstat-part-file",
        str(generated / "pars/simstat/partitions_ultem_insulators_Pb214_to_Po214.yaml"),
        "--detector-usabilities-file",
        str(generated / "pars/detector_usabilities.yaml"),
        "--simflow-config",
        str(config_path),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    hit.main()

    assert hit_file.exists(), "hit output file was not created"

    root_keys = lh5.ls(hit_file)
    assert "hit" in root_keys, f"'hit' group missing from {hit_file}; got {root_keys}"
    assert "tcm" in root_keys, f"'tcm' group missing from {hit_file}; got {root_keys}"

    # at least one detector table must be present under hit/
    det_tables = lh5.ls(hit_file, "hit/")
    assert len(det_tables) > 0, "no detector tables found under hit/"

    # verify the expected fields in the first detector table
    expected_fields = {
        "aoe",
        "aoe_raw",
        "drift_time_amax",
        "energy",
        "evtid",
        "is_single_site",
        "period",
        "psd_usability",
        "run",
        "t0",
        "usability",
    }
    first_det = det_tables[0]
    det_fields = set(lh5.ls(hit_file, first_det + "/"))
    # strip the leading path prefix from the returned names
    det_fields = {f.removeprefix(first_det + "/") for f in det_fields}
    missing = expected_fields - det_fields
    assert not missing, (
        f"fields {missing} missing from detector table {first_det}; got {det_fields}"
    )
