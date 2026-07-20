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
# way it's still possible to use Julia from a rule-specific conda env.
#
# Project.toml is an input so that editing it re-triggers instantiate; without it
# the input-less rule would be considered permanently up-to-date. Manifest.toml
# is gitignored (regenerated locally) so it must not be declared here.
rule _init_julia_env:
    localrule: True
    message:
        "Initializing Julia environment"
    input:
        "workflow/src/LegendSimflow.jl/Project.toml",
    output:
        config.paths.generated / ".julia-env-initialized",
    log:
        patterns.log_dirname(config) / "init-julia-env.log",
    shell:
        "julia --project=workflow/src/LegendSimflow.jl "
        "workflow/src/legendsimflow/scripts/init-julia-env.jl "
        "> {log} 2>&1 && touch {output}"


# the MAURINA gamma cascade data lives in a separate (currently private, ~1 GB) repository,
# Only simulations that opt in via the`maurina_gamma_cascades` simconfig field depend on this rule.
rule fetch_maurina_gamma_cascades:
    """Make the MAURINA gamma cascade data available on disk.

    Clones the `legend-exp/generated-maurina-output` repository into
    the writable, per-production ``paths.maurina_gamma_cascades`` directory. If the
    directory is already a git checkout, that revision is fetched and checked
    out.
    """
    localrule: True
    message:
        "Fetching MAURINA gamma cascade data"
    params:
        dest=config.paths.get("maurina_gamma_cascades", ""),
        rev=config.get("maurina_output_version", "main"),
        url="https://github.com/legend-exp/generated-maurina-output",
    output:
        touch(config.paths.generated / ".maurina-gamma-cascades-fetched"),
    log:
        patterns.log_dirname(config) / "fetch-maurina-gamma-cascades.log",
    shell:
        "{{ "
        "if [ -z '{params.dest}' ]; then "
        "  echo 'paths.maurina_gamma_cascades is not set in the Simflow config'; "
        "  exit 1; "
        "fi; "
        "if [ -d '{params.dest}/.git' ]; then "
        "  git -C '{params.dest}' fetch --depth 1 origin '{params.rev}' && "
        "  git -C '{params.dest}' checkout --detach FETCH_HEAD; "
        "else "
        "  git clone --depth 1 --single-branch --branch '{params.rev}' "
        "    '{params.url}' '{params.dest}'; "
        "fi; "
        "}} > {log} 2>&1"


# Memoize the on-demand fallback result to avoid re-running the expensive
# gen_list_of_all_hpges_valid_for_modeling() more than once per Snakemake
# invocation (touch executor / no YAML on disk).
HPGE_MODELING_CACHE: dict | None = None


def smk_load_hpge_cache() -> dict:
    """Load the modelable HPGe detectors, triggering the checkpoint if needed.

    Call this inside a Snakemake input function (lambda) or ``params:``
    lambda. Calling ``checkpoints.cache_modelable_hpges.get()`` triggers DAG
    re-evaluation: if the checkpoint has not run yet (or its output is stale),
    Snakemake will schedule it first and re-evaluate the calling rule's inputs
    once it completes.

    The on-disk cache lists every deployed HPGe with an ``is_modelable`` flag;
    this function returns only the modelable detectors (still keyed by name, per
    runid), so callers can expand the per-detector modeling rules directly.

    Falls back to computing the cache on-demand when running under the touch
    executor, which marks the checkpoint complete without executing its
    ``run:`` block (so the YAML may not exist on disk).
    """
    global HPGE_MODELING_CACHE

    outputs = checkpoints.cache_modelable_hpges.get().output
    is_modelable_path = Path(outputs.is_modelable)
    opv_path = Path(outputs.operational_voltage_in_V)
    if is_modelable_path.exists() and opv_path.exists():
        is_modelable = dbetto.utils.load_dict(is_modelable_path)
        opv = dbetto.utils.load_dict(opv_path)
    else:
        # touch executor marks the checkpoint complete without running it
        if HPGE_MODELING_CACHE is None:
            HPGE_MODELING_CACHE = aggregate.pivot_detinfo(
                aggregate.gen_list_of_all_hpges_valid_for_modeling(config)
            )
        is_modelable = HPGE_MODELING_CACHE["is_modelable"]
        opv = HPGE_MODELING_CACHE["operational_voltage_in_V"]

    return {
        runid: {
            hpge: {"operational_voltage_in_V": opv[runid][hpge]}
            for hpge, flag in dets.items()
            if flag
        }
        for runid, dets in is_modelable.items()
    }


checkpoint cache_modelable_hpges:
    """Cache the modeling status of the HPGe detectors.

    Computing which HPGe detectors are suitable for drift-time map and
    current-pulse model generation requires querying `legend-metadata` for
    every run in the Simflow. This can be slow, so the result is cached on disk
    under ``pars/detinfo/`` as one file per flag (``is_modelable.yaml`` and
    ``operational_voltage_in_V.yaml``), each a mapping ``runid -> detector ->
    value`` covering every deployed HPGe (both modelable and not). Downstream
    rules that need this information (e.g. `merge_hpge_drift_time_maps`) depend
    on this checkpoint so that the DAG can be re-evaluated after the cache is
    built.

    No wildcards are used.
    """
    localrule: True
    message:
        "Caching modelable HPGe detectors"
    params:
        runlist=config.runlist,
    output:
        is_modelable=patterns.detinfo_filename(config, "is_modelable"),
        operational_voltage_in_V=patterns.detinfo_filename(
            config, "operational_voltage_in_V"
        ),
    run:
        detinfo = aggregate.pivot_detinfo(
            aggregate.gen_list_of_all_hpges_valid_for_modeling(config)
        )
        dbetto.utils.write_dict(detinfo["is_modelable"], output.is_modelable)
        dbetto.utils.write_dict(
            detinfo["operational_voltage_in_V"], output.operational_voltage_in_V
        )


rule cache_detector_usabilities:
    """Cache detector usabilities.

    Querying the metadata for detector usability can be slow and constitute
    the bottleneck in post-processing (``opt`` and ``hit`` tiers). This rule
    caches, under ``pars/detinfo/``, one file per flag (``usability.yaml``,
    ``psd_usability.yaml``, ``crystal_metadata_usability.yaml``), each a mapping
    ``runid -> detector -> value``.
    """
    localrule: True
    message:
        "Caching detector usabilities"
    params:
        runlist=config.runlist,
    output:
        usability=patterns.detinfo_filename(config, "usability"),
        psd_usability=patterns.detinfo_filename(config, "psd_usability"),
        crystal_metadata_usability=patterns.detinfo_filename(
            config, "crystal_metadata_usability"
        ),
    run:
        detinfo = aggregate.pivot_detinfo(
            aggregate.gen_list_of_all_usabilities(config).to_dict()
        )
        dbetto.utils.write_dict(detinfo["usability"], output.usability)
        dbetto.utils.write_dict(detinfo["psd_usability"], output.psd_usability)
        dbetto.utils.write_dict(
            detinfo["crystal_metadata_usability"],
            output.crystal_metadata_usability,
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
        # deferred: avoid the expensive evaluation at Snakefile parse time
        lambda wc: aggregate.gen_list_of_all_plots(config),
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
