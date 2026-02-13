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

from legendsimflow import aggregate


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


# we use a dedicated dummy rule to initialize the Julia environment, in this
# way it's still possible to use Julia from a rule-specific conda env
rule _init_julia_env:
    localrule: True
    message:
        "Initializing Julia environment"
    output:
        config.paths.dtmaps / ".julia-env-initialized",
    shell:
        "cd workflow/src/legendsimflow/scripts && "
        "julia --project=. ./init-julia-env.jl && "
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
        config.paths.generated / "detector_usabilities.yaml",
    run:
        import dbetto

        dbetto.utils.write_dict(
            aggregate.gen_list_of_all_usabilities(config),
            output[0],
        )
