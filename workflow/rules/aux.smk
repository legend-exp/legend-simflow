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

from legendsimflow import aggregate, nersc


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
    the Simflow config field `make_tiers`.

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
    configuration.

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


rule make_sim_settings:
    """Suggest ``primaries_per_job`` and ``number_of_jobs`` from a benchmark run.

    Can be run with ``snakemake make_sim_settings``. Reads the benchmark log
    files produced by the benchmarking workflow and computes recommended
    values for ``primaries_per_job`` and ``number_of_jobs`` for each
    simulation, satisfying the following constraints:

    * No job runs longer than ``max_runtime_mins`` minutes.
    * The total simulated statistics exceed ``target_total_primaries``.
    * The total number of jobs is minimised.
    * Both values have at most ``n_sig_figs`` significant figures.

    The suggested settings are printed as a table and written to a YAML
    file (``{paths.generated}/sim_settings_suggestion.yaml``) that can be
    merged directly into the ``stp`` tier ``simconfig.yaml``.

    The constraints are controlled via the ``settings_opt`` section of the
    Simflow configuration:

    .. code-block:: yaml

        settings_opt:
          max_runtime_mins: 30          # default: 30
          target_total_primaries: 1e7   # default: 10_000_000
          n_sig_figs: 2                 # default: 2

    No wildcards are used.
    """
    localrule: True
    params:
        max_runtime_mins=config.get("settings_opt", {}).get("max_runtime_mins", 30),
        target_total_primaries=int(
            config.get("settings_opt", {}).get("target_total_primaries", 10_000_000)
        ),
        n_sig_figs=config.get("settings_opt", {}).get("n_sig_figs", 2),
    output:
        config.paths.generated / "sim_settings_suggestion.yaml",
    script:
        "../src/legendsimflow/scripts/make_sim_settings.py"


# we use a dedicated dummy rule to initialize the Julia environment, in this
# way it's still possible to use Julia from a rule-specific conda env
rule _init_julia_env:
    localrule: True
    message:
        "Initializing Julia environment"
    output:
        config.paths.generated / ".julia-env-initialized",
    shell:
        "julia --project=workflow/src/LegendSimflow.jl "
        "workflow/src/legendsimflow/scripts/init-julia-env.jl && "
        "touch {output}"


rule cache_detector_usabilities:
    """Cache detector usabilities.

    Querying the metadata for detector analysis usability can be slow and
    constitute the bottleneck in post-processing (`opt` and `hit` tiers).
    This rule caches the mapping `run -> detector -> usability` on disk.
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
