from legendsimflow import patterns, aggregate


rule gen_all_tier_evt:
    """Aggregate and produce all the `evt` tier files."""
    input:
        aggregate.gen_list_of_all_simid_outputs(config, tier="evt"),


rule build_tier_evt:
    """Produce an `evt` tier file.

    Event files re-organize the `hit` and `opt` tier data into a single,
    event-oriented table where each row correspond to an event.

    - a unified TCM is built from the `opt` and `hit` data. It is different
      from the `stp` tier TCM since it includes also the SiPM channels;
    - fields from lower tiers are restructured into events;
    - new event-level fields are computed and stored in the output file.

    Uses wildcards `simid` and `jobid`.
    """
    message:
        "Producing output file for job evt.{wildcards.simid}.{wildcards.jobid}"
    input:
        opt_file=patterns.output_simjob_filename(config, tier="opt"),
        hit_file=patterns.output_simjob_filename(config, tier="hit"),
    output:
        patterns.output_simjob_filename(config, tier="evt"),
    log:
        patterns.log_filename(config, SIMFLOW_CONTEXT.proctime, tier="evt"),
    benchmark:
        patterns.benchmark_filename(config, tier="evt")
    script:
        "../src/legendsimflow/scripts/tier/evt.py"
