from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest
from snakemake import api as smkapi

dummyprod = Path(__file__).parent / "dummyprod"
all_cores = os.cpu_count() or 1


def test_dag():
    output = smkapi.OutputSettings(verbose=False)

    # build workflow and DAG, execute with touch executor (no remage needed)
    with smkapi.SnakemakeApi(output) as api:
        wf_api = api.workflow(
            snakefile=dummyprod / "workflow/Snakefile",
            workdir=dummyprod,
            config_settings=smkapi.ConfigSettings(
                configfiles=(dummyprod / "simflow-config.yaml",)
            ),
            storage_settings=smkapi.StorageSettings(),
            resource_settings=smkapi.ResourceSettings(cores=all_cores),
        )
        dag = wf_api.dag()
        dag.execute_workflow(executor="touch")


def test_dag_simlist():
    # the explicit-simlist path must schedule the (PSD-gated) dtmap plots, which
    # the par tier contributes through process_simlist. The touch executor does
    # not create output files, so we check the scheduled jobs in the run log
    # instead (it also completes the modelable-HPGe checkpoint pulled in by the
    # requested hit tier, which a dry run would leave unresolved).
    output = smkapi.OutputSettings(verbose=False)

    with smkapi.SnakemakeApi(output) as api:
        wf_api = api.workflow(
            snakefile=dummyprod / "workflow/Snakefile",
            workdir=dummyprod,
            config_settings=smkapi.ConfigSettings(
                configfiles=(dummyprod / "simflow-config.yaml",),
                config={"simlist": "hit.pen_plates_Ra224_to_Pb208"},
            ),
            storage_settings=smkapi.StorageSettings(),
            resource_settings=smkapi.ResourceSettings(cores=all_cores),
        )
        dag = wf_api.dag()
        dag.execute_workflow(executor="touch")

    logs = sorted(
        (dummyprod / ".snakemake/log").glob("*.snakemake.log"),
        key=lambda p: p.stat().st_mtime,
    )
    assert "plot_hpge_drift_time_maps" in logs[-1].read_text(), (
        "dtmap plots not scheduled via the simlist path"
    )


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
