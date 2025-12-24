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

"""Rules to build the `stp` tier."""

from legendsimflow import aggregate, commands, utils, patterns
from legendsimflow import metadata as mutils


rule gen_all_tier_stp:
    """Build the entire `stp` tier."""
    input:
        aggregate.gen_list_of_all_plots_outputs(config, tier="stp"),
        aggregate.gen_list_of_all_simid_outputs(config, tier="stp"),


rule gen_geom_config:
    """Write a geometry configuration file for legend-pygeom-l200.

    Start from the template/default geometry configuration file and eventually
    add extra configuration options in case requested in `simconfig.yaml`
    through the `geom_config_extra` field.
    """
    message:
        "Generating geometry configuration for {wildcards.tier}.{wildcards.simid}"
    input:
        Path(config.paths.config) / "geom" / (config.experiment + "-geom-config.yaml"),
    params:
        # make this rule dependent on the actual simconfig block
        _simconfig_hash=lambda wc: mutils.smk_hash_simconfig(
            config, wc, "geom_config_extra"
        ),
    output:
        patterns.geom_config_filename(config),
    run:
        from dbetto import utils as dbetto_utils
        from legendsimflow import utils

        gconfig = dbetto_utils.load_dict(input[0])
        sconfig = mutils.get_simconfig(
            config, tier=wildcards.tier, simid=wildcards.simid
        )

        if "geom_config_extra" in sconfig:
            gconfig |= sconfig.geom_config_extra.to_dict()

        dbetto_utils.write_dict(gconfig, output[0])


rule build_geom_gdml:
    """Build a concrete geometry GDML file with {mod}`legend-pygeom-l200`.
    """
    message:
        "Building GDML geometry for {wildcards.tier}.{wildcards.simid}"
    input:
        rules.gen_geom_config.output,
    output:
        patterns.geom_gdml_filename(config),
    log:
        patterns.geom_log_filename(config, SIMFLOW_CONTEXT.proctime),
    shell:
        "LEGEND_METADATA={config.paths.metadata} "
        "legend-pygeom-l200 --verbose --config {input} -- {output} &> {log}"


def smk_remage_run(wildcards, input, output, threads):
    """Generate the remage command line for use in Snakemake rules."""
    return commands.remage_run(
        config,
        wildcards.simid,
        jobid=wildcards.jobid,
        tier="stp",
        geom=input.geom,
        output=output,
        procs=threads,
        macro_free=True,
    )


rule build_tier_stp:
    """Run a single simulation job for the `stp` tier.

    Uses wildcards `simid` and `jobid`.

    :::{note}
    The output remage file is declared as `protected` to avoid accidental
    deletions, since it typically takes a lot of resources to produce it.
    :::
    """
    message:
        "Producing output file for job stp.{wildcards.simid}.{wildcards.jobid}"
    input:
        verfile=lambda wc: patterns.vtx_filename_for_stp(config, wc.simid),
        geom=patterns.geom_gdml_filename(config, tier="stp"),
    params:
        cmd=smk_remage_run,
        # make this rule dependent on the actual simconfig block it is very
        # important here to ignore the simconfig fields that, if updated,
        # should not trigger re-creation of existing files. we ignore then
        # `number_of_jobs` and in the future we could also ignore
        # `primaries_per_job`, such that Snakemake will keep the existing stp
        # files with a different number of primaries on disk. Bonus: we ignore
        # `geom_config_extra` because that dependency is already tracked by
        # `input.geom`.
        _simconfig_hash=lambda wc: mutils.smk_hash_simconfig(
            config,
            wc,
            tier="stp",
            ignore=["geom_config_extra", "number_of_jobs"],
        ),
    output:
        # TODO: protected()
        patterns.output_simjob_filename(config, tier="stp"),
    log:
        patterns.log_filename(config, SIMFLOW_CONTEXT.proctime, tier="stp"),
    benchmark:
        patterns.benchmark_filename(config, tier="stp")
    threads: 1
    shell:
        # NOTE: since this can be a chain of commands, let's wrap it in {} to
        # make sure that all stderr/stdout is properly redirected to the log
        # file
        "{{ {params.cmd}; }} &> {log}"


rule plot_tier_stp_vertices:
    """Produces plots of the primary event vertices of tier `stp`.

    Only the first file of the simulation (i.e. job ID 0) is used. The rule
    is given a high priority to make sure that the plot is produced early. The
    maximum number of plotted events is set in the plotting script.

    Uses wildcard `simid`.
    """
    input:
        patterns.output_simjob_filename(config, tier="stp", jobid="0000"),
    output:
        patterns.plot_tier_stp_vertices_filename(config),
    priority: 100  # prioritize producing the needed input files over the others
    script:
        "../src/legendsimflow/scripts/plots/tier_stp_vertices.py"
