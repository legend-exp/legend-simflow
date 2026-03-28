from __future__ import annotations

from pathlib import Path

from snakemake import api as smkapi

dummyprod = Path(__file__).parent / "dummyprod"


def _run_workflow(extra_config=None):
    output = smkapi.OutputSettings(verbose=False)

    with smkapi.SnakemakeApi(output) as api:
        wf_api = api.workflow(
            snakefile=dummyprod / "workflow/Snakefile",
            workdir=dummyprod,
            config_settings=smkapi.ConfigSettings(
                configfiles=(dummyprod / "simflow-config.yaml",),
                config=extra_config or {},
            ),
            storage_settings=smkapi.StorageSettings(),
            resource_settings=smkapi.ResourceSettings(cores=1),
        )
        dag = wf_api.dag()
        dag.execute_workflow(executor="touch")


def test_run():
    _run_workflow()


def test_run_skip_psd():
    _run_workflow(extra_config={"options": {"skip_psd": True}})


def test_run_skip_lar():
    _run_workflow(extra_config={"options": {"skip_lar": True}})


def test_run_skip_psd_and_lar():
    _run_workflow(extra_config={"options": {"skip_psd": True, "skip_lar": True}})
