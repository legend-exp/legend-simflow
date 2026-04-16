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
        dag = wf_api.dag(smkapi.DAGSettings(forceall=True))
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
