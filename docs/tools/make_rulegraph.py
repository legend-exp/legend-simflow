"""Generate a Snakemake rule graph for the Simflow.

Runs in two phases:

1. Pre-computes the modelable-HPGe cache (normally produced by the
   ``cache_modelable_hpges`` checkpoint) so that post-checkpoint rules are
   visible to Snakemake's DAG builder.
2. Calls the Snakemake Python API to print the rule graph, with a monkey-patch
   for a Snakemake 9 bug where ``execution_settings`` is ``None`` in rulegraph
   mode after checkpoint resolution.

Usage (from the project root)::

    python docs/tools/make_rulegraph.py > docs/source/img/rulegraph.mmd

Or via the Makefile target::

    cd docs && make rulegraph
"""

from __future__ import annotations

import sys
from pathlib import Path

# snakemake.api must be imported before snakemake.dag to avoid a circular import
from snakemake import api as smkapi

# isort: split
import snakemake.dag as _smk_dag

ROOT = Path(__file__).resolve().parents[2]
DUMMYPROD = ROOT / "tests/dummyprod"
CONFIG = DUMMYPROD / "rulegraph-config.yaml"
SNAKEFILE = DUMMYPROD / "workflow/Snakefile"


def _guard_none_execution_settings(orig):
    def patched(self, job):
        if self.workflow.execution_settings is None:
            return False
        return orig(self, job)

    return patched


# workaround for snakemake bug: execution_settings is None in
# rulegraph mode after checkpoint resolution
_smk_dag.DAG.is_edit_notebook_job = _guard_none_execution_settings(
    _smk_dag.DAG.is_edit_notebook_job
)
_smk_dag.DAG.is_draft_notebook_job = _guard_none_execution_settings(
    _smk_dag.DAG.is_draft_notebook_job
)

sys.path.insert(0, str(ROOT / "workflow/src"))
import dbetto  # noqa: E402

from legendsimflow import aggregate, patterns  # noqa: E402
from legendsimflow.utils import init_simflow_context  # noqa: E402


def main() -> None:
    config = init_simflow_context(CONFIG).config

    # pre-write the modeling-status detinfo files so the checkpoint output
    # exists when the DAG is built
    detinfo = aggregate.pivot_detinfo(
        aggregate.gen_list_of_all_hpges_valid_for_modeling(config)
    )
    for flag in ("is_modelable", "operational_voltage_in_V"):
        path = patterns.detinfo_filename(config, flag)
        path.parent.mkdir(parents=True, exist_ok=True)
        dbetto.utils.write_dict(detinfo[flag], path)

    with smkapi.SnakemakeApi(smkapi.OutputSettings()) as api:
        dag_api = api.workflow(
            snakefile=SNAKEFILE,
            workdir=DUMMYPROD,
            config_settings=smkapi.ConfigSettings(configfiles=(CONFIG,)),
            storage_settings=smkapi.StorageSettings(),
            resource_settings=smkapi.ResourceSettings(cores=1),
        ).dag(dag_settings=smkapi.DAGSettings(print_dag_as="mermaid-js"))  # type: ignore[arg-type]
        dag_api.printrulegraph()


if __name__ == "__main__":
    main()
