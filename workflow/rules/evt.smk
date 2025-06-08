

rule gen_all_tier_evt:
    """Aggregate and produce all the evt tier files."""
    input:
        aggregate.gen_list_of_all_tier_evt_outputs(config),


rule make_tier_evt_config_file:
    """Generates configuration files for `build_tier_evt` based on metadata.
    Uses wildcard `runid`."""
    localrule: True
    input:
        # FIXME: need to list actual files, not the directory
        config["paths"]["metadata"],
    output:
        Path(config["paths"]["genconfig"]) / "{runid}-build_evt.json",
    script:
        "scripts/make_tier_evt_config_file.py"


rule make_run_partition_file:
    """Computes and stores on disk rules for partitioning the simulated event
    statistics according to data taking runs. Uses wildcard `simid`."""
    localrule: True
    input:
        hit_files=lambda wildcards: aggregate.gen_list_of_simid_outputs(
            config, tier="hit", simid=wildcards.simid
        ),
        runinfo=Path(config["paths"]["metadata"]) / "dataprod" / "runinfo.json",
    output:
        Path(config["paths"]["genconfig"]) / "{simid}-run_partition.json",
    params:
        ro_hit_files=lambda wildcards, input: utils.as_ro(config, input.hit_files),
    script:
        "scripts/make_run_partition_file.py"


rule build_tier_evt:
    """Produces an evt tier file."""
    message:
        "Producing output file for job evt.{wildcards.simid}.{wildcards.runid}"
    input:
        hit_files=lambda wildcards: aggregate.gen_list_of_simid_outputs(
            config, tier="hit", simid=wildcards.simid
        ),
        config_file=rules.make_tier_evt_config_file.output,
        run_part_file=rules.make_run_partition_file.output,
        hpge_db=Path(config["paths"]["metadata"])
        / "hardware/detectors/germanium/diodes",
    output:
        patterns.output_evt_filename(config),
    params:
        evt_window=lambda wildcards, input: tier_evt.smk_get_evt_window(
            wildcards, input
        ),
        hit_files_regex=utils.as_ro(
            config, patterns.output_simjob_regex(config, tier="hit")
        ),
    log:
        patterns.log_evtfile_path(config, proctime),
    benchmark:
        patterns.benchmark_evtfile_path(config)
    shell:
        patterns.run_command(config, "evt")
