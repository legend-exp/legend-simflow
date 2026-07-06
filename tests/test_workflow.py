from __future__ import annotations

import os
import shutil
from pathlib import Path

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


def _assert_psd_psl_in_hit(generated: Path) -> None:
    """Check the PSL-based detailed PSD output in the hit tier.

    A ``psd_psl`` sub-table must be produced for the modelable detector
    ``V05261B`` (it has a realistic pulse-shape library and drift-time map), and
    its A/E values must be at least partly finite, proving the detailed path
    actually ran rather than falling back to the all-NaN default.
    """
    hit_files = sorted((generated / "tier/hit").rglob("*.lh5"))
    assert hit_files, f"no hit output files found under {generated}/tier/hit"

    saw_psd_psl = False
    saw_finite_aoe = False
    for f in hit_files:
        if "hit/V05261B" not in lh5.ls(f, "hit/"):
            continue
        if "hit/V05261B/psd_psl" not in lh5.ls(f, "hit/V05261B/"):
            continue
        saw_psd_psl = True
        fields = {
            x.removeprefix("hit/V05261B/psd_psl/")
            for x in lh5.ls(f, "hit/V05261B/psd_psl/")
        }
        missing = _HIT_PSD_PSL_FIELDS - fields
        assert not missing, f"hit psd_psl fields {missing} missing in {f}; got {fields}"
        aoe = lh5.read_as("hit/V05261B/psd_psl/aoe", str(f), library="np")
        if np.any(np.isfinite(aoe)):
            saw_finite_aoe = True

    assert saw_psd_psl, "no hit file contains a hit/V05261B/psd_psl sub-table"
    assert saw_finite_aoe, (
        "hit/V05261B/psd_psl/aoe is all-NaN across all hit files; the PSL-based "
        "detailed PSD path did not compute real A/E values"
    )


def _assert_psd_psl_in_evt(generated: Path) -> None:
    """Check that the evt tier reads the detailed PSD back into geds/psd_psl."""
    evt_files = sorted((generated / "tier/evt").rglob("*.lh5"))
    assert evt_files, f"no evt output files found under {generated}/tier/evt"

    saw_psd_psl = False
    for f in evt_files:
        if "evt/geds" not in lh5.ls(f, "evt/"):
            continue
        if "evt/geds/psd_psl" not in lh5.ls(f, "evt/geds/"):
            continue
        saw_psd_psl = True
        fields = {
            x.removeprefix("evt/geds/psd_psl/") for x in lh5.ls(f, "evt/geds/psd_psl/")
        }
        missing = _EVT_PSD_PSL_FIELDS - fields
        assert not missing, f"evt psd_psl fields {missing} missing in {f}; got {fields}"

    assert saw_psd_psl, "no evt file contains an evt/geds/psd_psl sub-table"


# NOTE: the dry-run DAG-structure tests (DAG resolution, simlist scheduling,
# make_steps tier selection, PSD switches) live in test_dag.py.


@pytest.mark.needs_remage
@pytest.mark.skipif(shutil.which("remage") is None, reason="remage not installed")
def test_l1000_workflow():
    output = smkapi.OutputSettings(verbose=False)

    with smkapi.SnakemakeApi(output) as api:
        wf_api = api.workflow(
            snakefile=dummyprod / "workflow/Snakefile",
            workdir=dummyprod,
            config_settings=smkapi.ConfigSettings(
                configfiles=(dummyprod / "simflow-config-l1000.yaml",),
            ),
            storage_settings=smkapi.StorageSettings(),
            resource_settings=smkapi.ResourceSettings(cores=all_cores),
        )
        dag = wf_api.dag()
        dag.execute_workflow()

    # the l1000dsg01 hit settings enable simulate_psd_with_psl, so the detailed
    # (PSL-based) PSD must be produced in the hit tier and read back in the evt
    # tier
    generated = dummyprod / "generated-l1000"
    _assert_psd_psl_in_hit(generated)
    _assert_psd_psl_in_evt(generated)


@pytest.mark.needs_nersc
@pytest.mark.needs_remage
@pytest.mark.skipif(shutil.which("remage") is None, reason="remage not installed")
def test_l200_workflow():
    output = smkapi.OutputSettings(show_failed_logs=True)

    with smkapi.SnakemakeApi(output) as api:
        wf_api = api.workflow(
            snakefile=dummyprod / "workflow/Snakefile",
            workdir=dummyprod,
            config_settings=smkapi.ConfigSettings(
                configfiles=(dummyprod / "simflow-config-l200.yaml",),
            ),
            storage_settings=smkapi.StorageSettings(),
            resource_settings=smkapi.ResourceSettings(cores=all_cores),
        )
        dag = wf_api.dag()
        dag.execute_workflow(
            execution_settings=smkapi.ExecutionSettings(keep_going=True),
        )


@pytest.mark.needs_nersc
@pytest.mark.needs_remage
@pytest.mark.skipif(shutil.which("remage") is None, reason="remage not installed")
def test_l200cfg09_workflow():
    """Full vtx->pdf pipeline with the realistic (PSL-based) PSD on real data.

    Uses experiment ``l200cfg09`` with ``l200data`` at ref/v3.3.0 (NERSC-only)
    and a trimmed simconfig over two p16 ``ssc`` runs. PSD and the pulse-shape
    library are enabled, so this exercises the realistic-PSL-from-data path (the
    ``l1000dsg01`` test only builds the *ideal* PSL, having no ``l200data``).
    """
    output = smkapi.OutputSettings(show_failed_logs=True)

    with smkapi.SnakemakeApi(output) as api:
        wf_api = api.workflow(
            snakefile=dummyprod / "workflow/Snakefile",
            workdir=dummyprod,
            config_settings=smkapi.ConfigSettings(
                configfiles=(dummyprod / "simflow-config-l200cfg09.yaml",),
            ),
            storage_settings=smkapi.StorageSettings(),
            resource_settings=smkapi.ResourceSettings(cores=all_cores),
        )
        dag = wf_api.dag()
        dag.execute_workflow()

    # PSL is enabled with real l200data: the detailed (PSL-based) PSD must be
    # produced in the hit tier and read back in the evt tier
    generated = dummyprod / "generated-l200cfg09"
    _assert_psd_psl_in_hit(generated)
    _assert_psd_psl_in_evt(generated)
