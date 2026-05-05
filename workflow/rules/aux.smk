# Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from pathlib import Path

import dbetto
from legendsimflow import aggregate, nersc, patterns


def on_scratch_smk(path: str | Path):
    """Move `path` to the scratch directory if requested.

    This function is meant to be used in Snakemake rules to selectively make
    paths accessible from the scratch folder. If the NERSC scratch folder is
    disable in the Simflow config, returns `path` as is. Implemented via
    `snakemake-storage-plugin-fs`. Make sure to setup the local storage
    directory via Snakemake CLI or workflow profile.
    """
    if nersc.is_scratch_enabled(config) and not workflow.touch:
        return storage.fs(str(path))

    return path


rule print_stats:
    """Prints a table with summary runtime information for each `simid`.

    Can be run with `snakemake print_stats`. The listed tiers are taken from
    the Simflow config field `make_steps`.

    :::{note}
    The statistics refer to the total job wall time, as measured by Snakemake.
    :::

    No wildcards are used.
    """
    localrule: True
    script:
        "../src/legendsimflow/scripts/print_simflow_stats.py"


rule print_benchmark_stats:
    """Prints a table with summary runtime information of a benchmarking run.

    Can be run with `snakemake print_benchmark_stats`. This functionality is
    useful to tune the number of _remage_ primaries and jobs in the Simflow
    configuration. After printing the table, also writes an updated
    ``generated/benchmarks/generated-simconfig.yaml`` with suggested
    ``primaries_per_job`` and ``number_of_jobs`` values that can optionally
    be swapped in place of the source ``simconfig.yaml``.

    :::{note}
    The runtime and the simulation speed are extracted from the event
    simulation loop statistics reported by _remage_. These values do not
    account for other _remage_ steps like initialization or post-processing.
    :::

    No wildcards are used.
    """
    localrule: True
    script:
        "../src/legendsimflow/scripts/print_benchmark_stats.py"


# we use a dedicated dummy rule to initialize the Julia environment, in this
# way it's still possible to use Julia from a rule-specific conda env
rule _init_julia_env:
    localrule: True
    message:
        "Initializing Julia environment"
    output:
        config.paths.generated / ".julia-env-initialized",
    log:
        patterns.log_dirname(config) / "init-julia-env.log",
    shell:
        "julia --project=workflow/src/LegendSimflow.jl "
        "workflow/src/legendsimflow/scripts/init-julia-env.jl "
        "> {log} 2>&1 && touch {output}"


# Memoize the on-demand fallback result to avoid re-running the expensive
# gen_list_of_all_hpges_valid_for_modeling() more than once per Snakemake
# invocation (touch executor / no YAML on disk).
MODELABLE_HPGES: dict | None = None


def smk_load_hpge_cache() -> dict:
    """Load the modelable HPGe cache, triggering the checkpoint if needed.

    Call this inside a Snakemake input function (lambda) or ``params:``
    lambda. Calling ``checkpoints.cache_modelable_hpges.get()`` triggers DAG
    re-evaluation: if the checkpoint has not run yet (or its output is stale),
    Snakemake will schedule it first and re-evaluate the calling rule's inputs
    once it completes.

    Falls back to computing the cache on-demand when running under the touch
    executor, which marks the checkpoint complete without executing its
    ``run:`` block (so the YAML may not exist on disk).
    """
    global MODELABLE_HPGES

    yaml_path = Path(checkpoints.cache_modelable_hpges.get().output[0])
    if yaml_path.exists():
        return dbetto.utils.load_dict(yaml_path)
    # touch executor marks the checkpoint complete without running it
    if MODELABLE_HPGES is None:
        MODELABLE_HPGES = aggregate.gen_list_of_all_hpges_valid_for_modeling(config)
    return MODELABLE_HPGES


checkpoint cache_modelable_hpges:
    """Cache the list of HPGe detectors valid for modeling.

    Computing the list of HPGe detectors that are suitable for drift time map
    and current pulse model generation requires querying `legend-metadata` for
    every run in the Simflow. This can be slow, so the result is cached on disk
    as a YAML file mapping ``runid -> {hpge: voltage}``. Downstream rules that
    need this information (e.g. `merge_hpge_drift_time_maps`) depend on this
    checkpoint so that the DAG can be re-evaluated after the cache is built.

    No wildcards are used.
    """
    localrule: True
    message:
        "Caching modelable HPGe detectors"
    params:
        runlist=config.runlist,
    output:
        config.paths.pars / "modelable_hpge_detectors.yaml",
    run:
        aggregate.gen_list_of_all_hpges_valid_for_modeling(
            config, write_to_file=output[0]
        )


rule cache_detector_usabilities:
    """Cache detector usabilities.

    Querying the metadata for detector usability can be slow and constitute
    the bottleneck in post-processing (``opt`` and ``hit`` tiers). This rule
    caches the mapping ``run -> detector -> {usability, psd_usability}`` on
    disk.
    """
    localrule: True
    message:
        "Caching detector usabilities"
    params:
        runlist=config.runlist,
    output:
        config.paths.pars / "detector_usabilities.yaml",
    run:
        import dbetto

        dbetto.utils.write_dict(
            aggregate.gen_list_of_all_usabilities(config).to_dict(),
            output[0],
        )


rule archive_plots:
    """Archive all validation plots into a single tarball.

    Must be triggered manually with ``snakemake archive_plots`` — it is not
    part of the default ``all`` target. Collects all ``plots/`` subdirectories
    produced by the Simflow under the ``generated/`` directory and packs them
    into ``tarballs/<cycle>-plots.tar.xz``, preserving the directory tree
    structure.

    No wildcards are used.
    """
    localrule: True
    input:
        lambda wc: aggregate.gen_list_of_all_plots(config, cache=smk_load_hpge_cache()),
    output:
        patterns.plots_tarball_filename(config),
    run:
        from pathlib import Path
        from legendsimflow.archive import create_plots_tarball

        create_plots_tarball(
            generated_dir=config.paths.generated,
            output=Path(output[0]),
            prefix=Path(output[0]).name.removesuffix(".tar.xz"),
        )
