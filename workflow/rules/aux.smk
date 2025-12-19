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


rule print_stats:
    """Prints a table with summary runtime information for each `simid`.
    No wildcards are used."""
    localrule: True
    script:
        "../src/legendsimflow/scripts/print_simprod_stats.py"


rule print_benchmark_stats:
    """Prints a table with summary runtime information of a benchmarking run.
    No wildcards are used."""
    localrule: True
    script:
        "../src/legendsimflow/scripts/print_benchmark_stats.py"


# we use a dedicated dummy rule to initialize the Julia environment, in this
# way it's still possible to use Julia from a rule-specific conda env
rule _init_julia_env:
    message:
        "Initializing Julia environment"
    output:
        config.paths.dtmaps / ".julia-env-initialized",
    shell:
        "cd workflow/src/legendsimflow/scripts && "
        "julia --project=. ./init-julia-env.jl && "
        "touch {output}"
