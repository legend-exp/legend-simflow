from __future__ import annotations

import shutil
import sys
from pathlib import Path

import h5py
import lh5
import numpy as np
import pytest
import yaml

from legendsimflow.scripts.tier import hit

dummyprod = Path(__file__).parent.parent / "dummyprod"

_RUNIDS = ("l200-p03-r000-phy", "l200-p03-r001-phy")


@pytest.mark.needs_remage
def test_hit_script_cli(
    tmp_path,
    monkeypatch,
    legend_stp_path,
    legend_gdml_path,
    legend_dtmap_path,
    legend_currmod_paths,
    legend_hpge_obs_paths,
    legend_simstat_part_path,
    legend_detector_usabilities_path,
):
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

    dtmap_dir = pars_dir / "hpge/dtmaps"
    dtmap_dir.mkdir(parents=True)
    # r000: real dtmap for V05261B → PSD will be computed
    shutil.copy(legend_dtmap_path, dtmap_dir / f"{_RUNIDS[0]}-hpge-drift-time-maps.lh5")
    # r001: empty dtmap → dt_map = None → PSD will be NaN
    with h5py.File(dtmap_dir / f"{_RUNIDS[1]}-hpge-drift-time-maps.lh5", "w"):
        pass

    raw = yaml.safe_load((dummyprod / "simflow-config-l1000.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    raw["paths"]["pars"] = str(pars_dir)
    config_path = tmp_path / "simflow-config-l1000.yaml"
    config_path.write_text(yaml.safe_dump(raw))

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
            "--usability-file",
            str(legend_detector_usabilities_path / "usability.yaml"),
            "--psd-usability-file",
            str(legend_detector_usabilities_path / "psd_usability.yaml"),
            "--crystal-metadata-usability-file",
            str(legend_detector_usabilities_path / "crystal_metadata_usability.yaml"),
            "--simflow-config",
            str(config_path),
        ],
    )
    hit.main()

    assert hit_file.exists(), "hit output file was not created"

    root_keys = lh5.ls(hit_file)
    assert "hit" in root_keys, f"'hit' group missing from {hit_file}; got {root_keys}"
    assert "tcm" in root_keys, f"'tcm' group missing from {hit_file}; got {root_keys}"

    det_tables = lh5.ls(hit_file, "hit/")
    assert len(det_tables) > 0, "no detector tables found under hit/"

    expected_fields = {
        "psd",
        "energy",
        "evtid",
        "period",
        "psd_usability",
        "crystal_metadata_usability",
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

    # helper: read a single field from a detector table as a numpy array
    def _field(det: str, name: str) -> np.ndarray:
        return lh5.read_as(f"{det}/{name}", hit_file, library="np")

    # r000 has a dtmap for V05261B → finite PSD; r001 has none → NaN
    # both partitions land in the same hit/V05261B table, distinguishable by run
    assert "hit/V05261B" in det_tables, "V05261B not found in hit output"
    v_drift = _field("hit/V05261B/psd", "drift_time_amax")
    v_run = _field("hit/V05261B", "run")
    assert not np.all(np.isnan(v_drift[v_run == 0])), (
        "r000 V05261B drift_time_amax is all NaN despite having a dtmap"
    )
    assert np.all(np.isnan(v_drift[v_run == 1])), (
        "r001 V05261B drift_time_amax should be all NaN (no dtmap)"
    )

    # energy is always smeared — must be finite and positive
    energy = _field(first_det, "energy")
    assert np.all(np.isfinite(energy)), "energy contains non-finite values"
    assert np.all(energy >= 0), "energy contains negative values"

    # both runs are in period p03 → period == 3 throughout
    period_vals = _field(first_det, "period")
    assert np.all(period_vals == 3), (
        f"expected period=3 (p03), got unique={np.unique(period_vals)}"
    )

    # run values must come from the two configured runs (r000=0, r001=1)
    run_vals = _field(first_det, "run")
    assert np.all(np.isin(run_vals, [0, 1])), (
        f"unexpected run values: {np.unique(run_vals)}"
    )

    # usability and psd_usability must be valid encoded integers
    usability_vals = _field(first_det, "usability")
    assert np.all(np.isin(usability_vals, [0, 1, 2])), (
        f"invalid usability codes: {np.unique(usability_vals)}"
    )

    psd_usability_vals = _field(first_det, "psd_usability")
    assert np.all(np.isin(psd_usability_vals, [0, 1, 2])), (
        f"invalid psd_usability codes: {np.unique(psd_usability_vals)}"
    )

    crystal_metadata_usability_vals = _field(first_det, "crystal_metadata_usability")
    assert np.all(np.isin(crystal_metadata_usability_vals, [0, 1, 2])), (
        f"invalid crystal_metadata_usability codes: {np.unique(crystal_metadata_usability_vals)}"
    )
