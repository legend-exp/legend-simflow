
rule gen_all_tier_hit:
    """Aggregate and produce all the hit tier files."""
    input:
        aggregate.gen_list_of_all_simid_outputs(config, tier="hit"),


rule build_tier_hit:
    """Produces a hit tier file starting from a single stp tier file."""
    message:
        "Producing output file for job hit.{wildcards.simid}.{wildcards.jobid}"
    input:
        stp_file=patterns.input_simjob_filename(config, tier="hit"),
        optmap_lar=config["paths"]["optical_maps"]["lar"],
        optmap_pen=config["paths"]["optical_maps"]["pen"],
        optmap_fiber=config["paths"]["optical_maps"]["fiber"],
    output:
        patterns.output_simjob_filename(config, tier="hit"),
    params:
        ro_stp_file=lambda wildcards, input: utils.as_ro(config, input.stp_file),
    log:
        patterns.log_file_path(config, proctime, tier="hit"),
    benchmark:
        patterns.benchmark_file_path(config, tier="hit")
    shell:
        patterns.run_command(config, "hit")
