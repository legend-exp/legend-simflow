from legendsimflow import patterns, aggregate


rule gen_all_tier_evt:
    """Aggregate and produce all the `evt` tier files."""
    input:
        aggregate.gen_list_of_all_simid_outputs(config, tier="evt"),
        aggregate.gen_list_of_all_plots_outputs(config, tier="evt"),


rule build_tier_evt:
    """Produce an `evt` tier file.

    Event files re-organize the `hit` and `opt` tier data into a single,
    event-oriented table where each row correspond to an event.

    - a unified TCM is built from the `opt` and `hit` data. It is different
      from the `stp` tier TCM since it includes also the SiPM channels;
    - each chunk of the unified TCM is partitioned according to the livetime
      span of each run (see the `make_simstat_partition_file` rule);
    - fields from lower tiers are restructured into events;
    - new event-level fields are computed and stored in the output file;
    - optionally, random-coincidence (RC) SiPM data from real evt files is
      added as `spms/rc_energy` and `spms/rc_time` (controlled by
      ``add_random_coincidences`` in ``tier/evt/{experiment}/settings.yaml``).

    Uses wildcards `simid` and `jobid`.
    """
    message:
        "Producing output file for job evt.{wildcards.simid}.{wildcards.jobid}"
    input:
        stp_file=patterns.output_simjob_filename(config, tier="stp"),
        opt_file=patterns.output_simjob_filename(config, tier="opt"),
        hit_file=patterns.output_simjob_filename(config, tier="hit"),
        simstat_part_file=patterns.simstat_part_filename(config),
        detector_usabilities=rules.cache_detector_usabilities.output,
    params:
        add_random_coincidences=lambda wc: config.metadata.simprod.config.tier.evt[
            config.experiment
        ].settings.add_random_coincidences,
    output:
        patterns.output_simjob_filename(config, tier="evt"),
    log:
        patterns.log_filename(config, tier="evt"),
    benchmark:
        patterns.benchmark_filename(config, tier="evt")
    script:
        "../src/legendsimflow/scripts/tier/evt.py"
