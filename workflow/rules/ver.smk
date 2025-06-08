simconfigs = aggregate.collect_simconfigs(config, ["ver"])
for tier, simid, n_macros in simconfigs:

    rule:
        """Generates all needed simulation macros. No wildcards are used."""  # ({n_macros}) for {simid} in tier {tier}"""
        localrule: True
        input:
            **patterns.macro_gen_inputs(config, tier, simid),
        output:
            patterns.input_simid_filenames(config, n_macros, tier=tier, simid=simid),
        params:
            tier=tier,
            simid=simid,
        threads: 1
        message:
            f"Generating macros for {tier}.{simid}"
        script:
            "scripts/generate_macros.py"

    utils.set_last_rule_name(workflow, f"gen_macros_{simid}-tier_{tier}")

    rule:
        """Produces plots for the primary event vertices of simid {simid} in tier {tier}"""
        input:
            aggregate.gen_list_of_simid_outputs(config, tier, simid, max_files=5),
        output:
            Path(patterns.plots_file_path(config, tier=tier, simid=simid))
            / f"mage-event-vertices-tier_{tier}.png",
        priority: 100  # prioritize producing the needed input files over the others
        shell:
            (
                " ".join(config["execenv"])
                + " python "
                + workflow.source_path("../scripts/plot_mage_vertices.py")
                + " -b -o {output} {input}"
            )

    utils.set_last_rule_name(workflow, f"plot_prim_vert_{simid}-tier_{tier}")


rule build_tier_ver:
    """Run a single simulation job for the ver tier.
    Uses wildcards `simid` and `jobid`.

    Warning
    -------
    The macro file is marked as "ancient" as a workaround to the fact that
    it might have been re-generated (i.e. it effectively has a more recent
    creation time) but with the same content as before (i.e. there is no need
    to re-run the simulation). If the macro content is updated, users will need
    to manually remove the output simulation files or force execution.
    """
    message:
        "Producing output file for job ver.{wildcards.simid}.{wildcards.jobid}"
    input:
        macro=ancient(patterns.input_simjob_filename(config, tier="ver")),
    output:
        protected(patterns.output_simjob_filename(config, tier="ver")),
    log:
        patterns.log_file_path(config, proctime, tier="ver"),
    benchmark:
        patterns.benchmark_file_path(config, tier="ver")
    shell:
        patterns.run_command(config, "ver")
