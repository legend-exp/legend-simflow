from __future__ import annotations

import os
import shutil
from pathlib import Path

import awkward as ak
import dbetto
import lh5
import numpy as np
import pytest
from snakemake import api as smkapi

dummyprod = Path(__file__).parent / "dummyprod"
all_cores = os.cpu_count() or 1

# fields written by the PSL-based (detailed) PSD path
_HIT_PSD_PSL_FIELDS = {
    "drift_time_amax",
    "aoe_raw",
    "aoe_corr",
    "aoe",
    "is_single_site",
    "is_bb_like",
    "is_high_aoe",
}
_EVT_PSD_PSL_FIELDS = {
    "aoe",
    "has_aoe",
    "aoe_corr",
    "drift_time_amax",
    "is_single_site",
    "is_bb_like",
    "is_high_aoe",
}


def _assert_psd_psl_in_hit(generated: Path, det: str = "V00001A") -> None:
    """Check the PSL-based detailed PSD output in the hit tier.

    A ``psd/pulse_lib`` sub-table must be produced for the modelable detector
    *det* (it has a realistic pulse-shape library and drift-time map), and its
    A/E values must be at least partly finite, proving the detailed path
    actually ran rather than falling back to the all-NaN default.
    """
    hit_files = sorted((generated / "tier/hit").rglob("*.lh5"))
    assert hit_files, f"no hit output files found under {generated}/tier/hit"

    saw_psd_psl = False
    saw_finite_aoe = False
    for f in hit_files:
        if f"hit/{det}" not in lh5.ls(f, "hit/"):
            continue
        if f"hit/{det}/psd/pulse_lib" not in lh5.ls(f, f"hit/{det}/psd/"):
            continue
        saw_psd_psl = True
        fields = {
            x.removeprefix(f"hit/{det}/psd/pulse_lib/")
            for x in lh5.ls(f, f"hit/{det}/psd/pulse_lib/")
        }
        missing = _HIT_PSD_PSL_FIELDS - fields
        assert not missing, f"hit psd_psl fields {missing} missing in {f}; got {fields}"
        aoe = lh5.read_as(f"hit/{det}/psd/pulse_lib/aoe", str(f), library="np")
        if np.any(np.isfinite(aoe)):
            saw_finite_aoe = True

    assert saw_psd_psl, f"no hit file contains a hit/{det}/psd/pulse_lib sub-table"
    assert saw_finite_aoe, (
        f"hit/{det}/psd/pulse_lib/aoe is all-NaN across all hit files; the "
        "PSL-based detailed PSD path did not compute real A/E values"
    )


def _assert_psd_psl_in_evt(generated: Path) -> None:
    """Check that the evt tier reads the detailed PSD into geds/psd/pulse_lib."""
    evt_files = sorted((generated / "tier/evt").rglob("*.lh5"))
    assert evt_files, f"no evt output files found under {generated}/tier/evt"

    saw_psd_psl = False
    for f in evt_files:
        if "evt/geds" not in lh5.ls(f, "evt/"):
            continue
        if "evt/geds/psd/pulse_lib" not in lh5.ls(f, "evt/geds/psd/"):
            continue
        saw_psd_psl = True
        fields = {
            x.removeprefix("evt/geds/psd/pulse_lib/")
            for x in lh5.ls(f, "evt/geds/psd/pulse_lib/")
        }
        missing = _EVT_PSD_PSL_FIELDS - fields
        assert not missing, f"evt psd_psl fields {missing} missing in {f}; got {fields}"

    assert saw_psd_psl, "no evt file contains an evt/geds/psd/pulse_lib sub-table"


def _run_workflow(configfile: Path, keep_going: bool = False) -> None:
    """Run the dummy production configured in *configfile* on all cores."""
    output = smkapi.OutputSettings(verbose=False, show_failed_logs=True)

    with smkapi.SnakemakeApi(output) as api:
        wf_api = api.workflow(
            snakefile=dummyprod / "workflow/Snakefile",
            workdir=dummyprod,
            config_settings=smkapi.ConfigSettings(configfiles=(configfile,)),
            storage_settings=smkapi.StorageSettings(),
            resource_settings=smkapi.ResourceSettings(cores=all_cores),
        )
        dag = wf_api.dag()
        dag.execute_workflow(
            execution_settings=smkapi.ExecutionSettings(keep_going=keep_going),
        )


# NOTE: the dry-run DAG-structure tests (DAG resolution, simlist scheduling,
# make_steps tier selection, PSD switches) live in test_dag.py.


@pytest.mark.needs_remage
@pytest.mark.skipif(shutil.which("remage") is None, reason="remage not installed")
def test_l200_workflow():
    _run_workflow(dummyprod / "simflow-config-l200.yaml")

    # the l200cfg01 hit settings enable simulate_psd_with_psl, so the detailed
    # (PSL-based) PSD must be produced in the hit tier and read back in the evt
    # tier
    generated = dummyprod / "generated-l200"
    _assert_psd_psl_in_hit(generated)
    _assert_psd_psl_in_evt(generated)

    # the SSD simulations are tuned to the depletion voltages (in V) of the
    # mock impurity model in paths.hpge_impurity_model
    ssd_info = dbetto.utils.load_dict(generated / "pars/detinfo/hpge_ssd_modeling.yaml")
    tuned = {"V00001A": 3000, "V00001B": 3050}
    for dets in ssd_info.values():
        for det, info in dets.items():
            assert info["measured_depletion_voltage_in_V"] == tuned[det]


@pytest.mark.needs_remage
@pytest.mark.skipif(shutil.which("remage") is None, reason="remage not installed")
def test_l1000_workflow():
    _run_workflow(dummyprod / "simflow-config-l1000.yaml")

    # same PSD settings as l200cfg01: the detailed (PSL-based) PSD must be
    # produced in the hit tier and read back in the evt tier
    generated = dummyprod / "generated-l1000"
    _assert_psd_psl_in_hit(generated, det="V00101Z")
    _assert_psd_psl_in_evt(generated)

    # the HPGe detectors come from the legend-pygeom-l1000 configuration (one
    # string of two), their records are derived from the default record
    ssd_info = dbetto.utils.load_dict(generated / "pars/detinfo/hpge_ssd_modeling.yaml")
    for dets in ssd_info.values():
        assert set(dets) == {"V00101Z", "V00102Z"}

    # the detector identifiers in the evt tier are the ones pylegendmeta
    # derives for string 1, positions 1 and 2
    evt_files = sorted((generated / "tier/evt").rglob("*.lh5"))
    assert evt_files, f"no evt output files found under {generated}/tier/evt"
    for f in evt_files:
        rawid = lh5.read_as("evt/geds/rawid", str(f), library="ak")
        assert set(np.unique(ak.flatten(rawid).to_numpy())) <= {1001001, 1001002}

    # the optics are simulated with the combined optical map
    opt_files = sorted((generated / "tier/opt").rglob("*.lh5"))
    assert opt_files, f"no opt output files found under {generated}/tier/opt"
    for f in opt_files:
        assert "hit/spms" in lh5.ls(f, "hit/")


@pytest.mark.needs_nersc
@pytest.mark.needs_remage
@pytest.mark.skipif(shutil.which("remage") is None, reason="remage not installed")
def test_l200_nersc_workflow():
    _run_workflow(dummyprod / "simflow-config-l200-nersc.yaml", keep_going=True)
