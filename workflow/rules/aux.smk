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
from legendsimflow.exceptions import SimflowConfigError
from legendsimflow.metadata import (
    get_tier_settings,
    get_par_settings,
    simulates_hpge_psd,
)

# the HPGe modeling cache is only needed by the PSD-gated par outputs. aux.smk
# is always included, so the flag is available to the other rule modules too
_simulate_psd = simulates_hpge_psd(config)
_tune_impurity = get_tier_settings(config, "hit").tune_hpge_impurities_on_data
_impurity_settings = get_par_settings(config, "impurityscan")
# a fit window starting below the scan cut would compare data with an
# incomplete simulated distribution
if _tune_impurity and (
    _impurity_settings.drift_time_scan.min_energy_in_keV
    > _impurity_settings.impurity_fit.energy_range_in_keV[0]
):
    raise SimflowConfigError(
        "drift_time_scan.min_energy_in_keV is above the lower edge of "
        "impurity_fit.energy_range_in_keV",
        "pars.impurityscan.settings",
    )

if _tune_impurity and not _simulate_psd:
    raise SimflowConfigError(
        "tune_hpge_impurities_on_data requires simulate_psd", "tier.hit.settings"
    )
if _tune_impurity and config.paths.get("l200data", None) is None:
    raise SimflowConfigError(
        "tune_hpge_impurities_on_data requires the LEGEND-200 data (paths.l200data)",
        "simflow-config.paths",
    )


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


# loaded once per Snakemake invocation, see smk_load_hpge_cache()
HPGE_MODELING_CACHE: dict | None = None


def smk_load_hpge_cache() -> dict:
    """Load the modelable HPGe detectors, triggering the checkpoint if needed.

    Call this inside a Snakemake input function (lambda) or ``params:``
    lambda. Calling ``checkpoints.cache_modelable_hpges.get()`` triggers DAG
    re-evaluation: if the checkpoint has not run yet (or its output is stale),
    Snakemake will schedule it first and re-evaluate the calling rule's inputs
    once it completes. Once ``get()`` returns, the checkpoint output does not
    change for the rest of the Snakemake invocation, so it is loaded only once.

    The on-disk cache lists every deployed HPGe with an ``is_modelable`` flag;
    this function returns only the modelable detectors (still keyed by name, per
    runid), so callers can expand the per-detector modeling rules directly.

    Falls back to computing the cache on-demand when running under the touch
    executor, which marks the checkpoint complete without executing its
    ``run:`` block (so the YAML may not exist on disk).
    """
    global HPGE_MODELING_CACHE

    outputs = checkpoints.cache_modelable_hpges.get().output
    if HPGE_MODELING_CACHE is not None:
        return HPGE_MODELING_CACHE

    is_modelable_path = Path(outputs.is_modelable)
    opv_path = Path(outputs.operational_voltage_in_V)
    if is_modelable_path.exists() and opv_path.exists():
        is_modelable = dbetto.utils.load_dict(is_modelable_path)
        opv = dbetto.utils.load_dict(opv_path)
    else:
        # touch executor marks the checkpoint complete without running it
        detinfo = aggregate.pivot_detinfo(
            aggregate.gen_list_of_all_hpges_valid_for_modeling(config)
        )
        is_modelable = detinfo["is_modelable"]
        opv = detinfo["operational_voltage_in_V"]

    HPGE_MODELING_CACHE = {
        runid: {
            hpge: {"operational_voltage_in_V": opv[runid][hpge]}
            for hpge, flag in dets.items()
            if flag
        }
        for runid, dets in is_modelable.items()
    }
    return HPGE_MODELING_CACHE


def smk_wait_for_hpge_cache(wildcards) -> list:
    """Make a job wait for the `cache_modelable_hpges` checkpoint.

    For rules whose ``params`` read the cache. Snakemake requires the job input
    to depend on the checkpoint, and the cache files themselves as input would
    rerun the job whenever any detector changes.
    """
    smk_load_hpge_cache()
    return []


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

    The cache is rebuilt when the runs, the modeling settings or the tuned HPGe
    impurity model (``paths.hpge_impurity_model``) of the Simflow change, but
    not when the metadata changes: the metadata files are not
    inputs, since a change in any of them, even for other detectors or runs,
    would rebuild the cache. After a metadata update, force the rebuild with
    ``--forcerun cache_modelable_hpges``.

    No wildcards are used.
    """
    localrule: True
    message:
        "Caching modelable HPGe detectors"
    input:
        impurity=config.paths.get("hpge_impurity_model", []),
    params:
        runids=lambda wc: sorted(aggregate.gen_list_of_all_runids(config)),
        tune_hpge_impurities_on_data=_tune_impurity,
        min_voltage_above_depletion_in_V=get_par_settings(
            config, "modeling"
        ).min_voltage_above_depletion_in_V,
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
    ``rawid.yaml``, ``psd_usability.yaml``,
    ``crystal_metadata_usability.yaml``), each a mapping
    ``runid -> detector -> value``.

    The cache is rebuilt when the runs of the Simflow change, but not when the
    metadata changes (see `cache_modelable_hpges`). After a metadata update,
    force the rebuild with ``--forcerun cache_detector_usabilities``.
    """
    localrule: True
    message:
        "Caching detector usabilities"
    params:
        runids=lambda wc: sorted(aggregate.gen_list_of_all_runids(config)),
        random_coincidence_runid=get_tier_settings(
            config, "evt"
        ).random_coincidence_runid,
    output:
        usability=patterns.detinfo_filename(config, "usability"),
        rawid=patterns.detinfo_filename(config, "rawid"),
        psd_usability=patterns.detinfo_filename(config, "psd_usability"),
        crystal_metadata_usability=patterns.detinfo_filename(
            config, "crystal_metadata_usability"
        ),
    run:
        detinfo = aggregate.pivot_detinfo(
            aggregate.gen_list_of_all_usabilities(config).to_dict()
        )
        dbetto.utils.write_dict(detinfo.get("usability", {}), output.usability)
        dbetto.utils.write_dict(detinfo.get("rawid", {}), output.rawid)
        dbetto.utils.write_dict(detinfo.get("psd_usability", {}), output.psd_usability)
        dbetto.utils.write_dict(
            detinfo.get("crystal_metadata_usability", {}),
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
        # don't trigger the checkpoint if no par-tier plot is part of the build
        lambda wc: aggregate.gen_list_of_all_plots(
            config,
            cache=(
                smk_load_hpge_cache()
                if _simulate_psd and "par" in config.make_steps
                else None
            ),
        ),
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
