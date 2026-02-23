import snakemake.io

from legendsimflow import metadata as mutils
from legendsimflow import patterns


def smk_build_tier_vtx_shell_cmd(wildcards, input, output):
    """Generate the `build_tier_vtx` Snakemake rule command line."""
    n_events = mutils.get_simconfig(config, "stp", wildcards.simid, "primaries_per_job")

    # expand the wildcards contained in the vtx-generation command provided by the user
    return snakemake.io.expand(
        mutils.get_vtx_simconfig(config, wildcards.simid).command,
        GDML_FILE=input.geom,
        OUTPUT_FILE=output,
        N_EVENTS=n_events,
    )


rule build_tier_vtx:
    """Run a single simulation job for the vtx tier.

    Uses wildcards `simid` and `jobid`.
    """
    message:
        "Producing output file for job vtx.{wildcards.simid}.{wildcards.jobid}"
    input:
        geom=on_scratch_smk(patterns.geom_gdml_filename(config, tier="stp")),
    params:
        # NOTE: a change in the simconfig command will be detected by Snakemake
        cmd=smk_build_tier_vtx_shell_cmd,
    output:
        patterns.output_simjob_filename(config, tier="vtx"),
    log:
        patterns.log_filename(config, tier="vtx"),
    benchmark:
        patterns.benchmark_filename(config, tier="vtx")
    shell:
        "{params.cmd} &> {log}"
