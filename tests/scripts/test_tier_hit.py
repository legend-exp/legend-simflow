from __future__ import annotations

import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest
import yaml
from lgdo import lh5

from legendsimflow.scripts.tier import hit

dummyprod = Path(__file__).parent.parent / "dummyprod"

_RUNIDS = ("l200-p03-r000-phy", "l200-p03-r001-phy")


@pytest.mark.needs_remage
def test_hit_script_cli(
    tmp_path,
    monkeypatch,
    legend_stp_path,
    legend_gdml_path,
    legend_currmod_paths,
    legend_hpge_obs_paths,
    legend_simstat_part_path,
    legend_detector_usabilities_path,
):
    # ------------------------------------------------------------------ pars
    pars_dir = tmp_path / "pars"

    for runid in _RUNIDS:
        obs = legend_hpge_obs_paths[runid]
        for subdir, src, dest_name in (
            ("hpge/eresmod", obs["eresmod"], f"{runid}-model.yaml"),
            ("hpge/aoeresmod", obs["aoeresmod"], f"{runid}-model.yaml"),
            ("hpge/psdcuts", obs["psdcuts"], f"{runid}-psd-cuts.yaml"),
            ("hpge/currmod", legend_currmod_paths[runid], f"{runid}-model.yaml"),
        ):
            dest_dir = pars_dir / subdir
            dest_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy(src, dest_dir / dest_name)

    # empty dtmap files — file must exist; no valid maps → dt_map = None (no crash)
    dtmap_dir = pars_dir / "hpge/dtmaps"
    dtmap_dir.mkdir(parents=True)
    for runid in _RUNIDS:
        with h5py.File(dtmap_dir / f"{runid}-hpge-drift-time-maps.lh5", "w"):
            pass

    # --------------------------------------------------------- simflow config
    raw = yaml.safe_load((dummyprod / "simflow-config-l1000.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    raw["paths"]["pars"] = str(pars_dir)
    config_path = tmp_path / "simflow-config-l1000.yaml"
    config_path.write_text(yaml.safe_dump(raw))

    # ---------------------------------------------------------------- run hit
    hit_file = tmp_path / "hit.lh5"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "hit",
            "--stp-file",
            str(legend_stp_path),
            "--jobid",
            "0000",
            "--hit-file",
            str(hit_file),
            "--geom-file",
            str(legend_gdml_path),
            "--simstat-part-file",
            str(legend_simstat_part_path),
            "--detector-usabilities-file",
            str(legend_detector_usabilities_path),
            "--simflow-config",
            str(config_path),
        ],
    )
    hit.main()

    # -------------------------------------------------------------- assertions
    assert hit_file.exists(), "hit output file was not created"

    root_keys = lh5.ls(hit_file)
    assert "hit" in root_keys, f"'hit' group missing from {hit_file}; got {root_keys}"
    assert "tcm" in root_keys, f"'tcm' group missing from {hit_file}; got {root_keys}"

    det_tables = lh5.ls(hit_file, "hit/")
    assert len(det_tables) > 0, "no detector tables found under hit/"

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
    det_fields = {
        f.removeprefix(first_det + "/") for f in lh5.ls(hit_file, first_det + "/")
    }
    missing = expected_fields - det_fields
    assert not missing, (
        f"fields {missing} missing from detector table {first_det}; got {det_fields}"
    )

    # helper: read a single field from the first detector table as a numpy array
    def _field(name: str) -> np.ndarray:
        return lh5.read_as(f"{first_det}/{name}", hit_file, library="np")

    # empty dtmap files → dt_map is None → PSD observables are all NaN
    assert np.all(np.isnan(_field("drift_time_amax"))), (
        "drift_time_amax should be all NaN (no dtmap)"
    )
    assert np.all(np.isnan(_field("aoe_raw"))), "aoe_raw should be all NaN (no dtmap)"
    assert np.all(np.isnan(_field("aoe"))), "aoe should be all NaN (no dtmap)"
    assert not np.any(_field("is_single_site")), (
        "is_single_site should be all False (no dtmap)"
    )

    # energy is always smeared — must be finite and positive
    energy = _field("energy")
    assert np.all(np.isfinite(energy)), "energy contains non-finite values"
    assert np.all(energy >= 0), "energy contains negative values"

    # both runs are in period p03 → period == 3 throughout
    period_vals = _field("period")
    assert np.all(period_vals == 3), (
        f"expected period=3 (p03), got unique={np.unique(period_vals)}"
    )

    # run values must come from the two configured runs (r000=0, r001=1)
    run_vals = _field("run")
    assert np.all(np.isin(run_vals, [0, 1])), (
        f"unexpected run values: {np.unique(run_vals)}"
    )

    # usability and psd_usability must be valid encoded integers
    usability_vals = _field("usability")
    assert np.all(np.isin(usability_vals, [0, 1, 2])), (
        f"invalid usability codes: {np.unique(usability_vals)}"
    )

    psd_usability_vals = _field("psd_usability")
    assert np.all(np.isin(psd_usability_vals, [0, 1, 2])), (
        f"invalid psd_usability codes: {np.unique(psd_usability_vals)}"
    )
