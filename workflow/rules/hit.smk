from legendsimflow import aggregate, patterns


rule gen_all_tier_hit:
    """Aggregate and produce all the hit tier files."""
    input:
        aggregate.gen_list_of_all_simid_outputs(config, tier="hit"),
        lambda wc: aggregate.gen_list_of_all_plots_outputs(
            config, tier="hit", cache=smk_load_hpge_cache()
        ),


# NOTE: we don't rely on rules from other tiers here (e.g.
# rules.build_tiers_stp.output) because we want to support making only the hit
# tier via the config.make_steps option
rule build_tier_hit:
    """Produce a `hit` tier file starting from a single `stp` tier file.

    This rule implements the post-processing of the `stp` tier HPGe data in
    chunks, in the following steps:

    - each chunk is partitioned according to the livetime span of each run
      (see the `make_simstat_partition_file` rule). For each partition:
    - the detector usability and PSD usability are retrieved from
      `legend-metadata` and stored in the output;
    - the active volume model is applied based on information from `legend-metadata`;
    - A/E is simulated based on current signal templates extracted from
      LEGEND-200 data;
    - energy is smeared according to the measured energy resolution (extracted
      from the data production parameters database);
    - a new time-coincidence map (TCM) across the processed detectors is
      created and stored in the output file.

    The `stp` data format is preserved: detector tables are stored separately
    in the output file below `/hit/{detector_name}`.

    Uses wildcards `simid` and `jobid`.
    """
    message:
        "Producing output file for job hit.{wildcards.simid}.{wildcards.jobid}"
    input:
        geom=patterns.geom_gdml_filename(config, tier="stp"),
        stp_file=patterns.output_simjob_filename(config, tier="stp"),
        # NOTE: we pass here the full list of maps/current models, but likely
        # not all of them will be used. room for improvement
        hpge_dtmaps=lambda wc: aggregate.gen_list_of_merged_dtmaps(config, wc.simid),
        hpge_currmods=lambda wc: aggregate.gen_list_of_merged_currmods(config, wc.simid),
        hpge_eresmods=lambda wc: aggregate.gen_list_of_eresmods(config, wc.simid),
        hpge_aoeresmods=lambda wc: aggregate.gen_list_of_aoeresmods(config, wc.simid),
        hpge_psdcuts=lambda wc: aggregate.gen_list_of_psdcuts(config, wc.simid),
        # NOTE: technically this rule only depends on one block in the
        # partitioning file, but in practice the full file will always change
        simstat_part_file=patterns.simstat_part_filename(config),
        detector_usabilities=rules.cache_detector_usabilities.output,
    output:
        patterns.output_simjob_filename(config, tier="hit"),
    log:
        patterns.log_filename(config, tier="hit"),
    benchmark:
        patterns.benchmark_filename(config, tier="hit")
    script:
        "../src/legendsimflow/scripts/tier/hit.py"


rule plot_tier_hit_observables:
    """Produce validation plots of observable distributions from the `hit` tier.

    Generates diagnostic plots from all `hit` output files for the given
    `simid`.

    Uses wildcard `simid`.
    """
    message:
        "Producing control plots for job hit.{wildcards.simid}"
    input:
        lambda wc: aggregate.gen_list_of_simid_outputs(
            config, tier="hit", simid=wc.simid
        ),
    output:
        patterns.plot_tier_hit_observables_filename(config),
    script:
        "../src/legendsimflow/scripts/plots/tier_hit_observables.py"
