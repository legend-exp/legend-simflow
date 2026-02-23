from legendsimflow import patterns, aggregate


rule gen_all_tier_cvt:
    """Aggregate and produce all the `cvt` tier files."""
    input:
        aggregate.gen_list_of_all_tier_cvt_outputs(config),
        aggregate.gen_list_of_all_plots_outputs(config, tier="cvt"),


rule build_tier_cvt:
    """Produce a `cvt` tier file.

    `cvt` stands for "concatenated `evt` tier". `evt` files for each simulation
    job are concatenated/aggregated into a single file.

    Uses wildcards `simid`.
    """
    message:
        "Producing output file for job cvt.{wildcards.simid}"
    input:
        lambda wildcards: aggregate.gen_list_of_simid_outputs(
            config, tier="evt", simid=wildcards.simid
        ),
    output:
        patterns.output_tier_cvt_filename(config),
    log:
        patterns.log_tier_cvt_filename(config),
    benchmark:
        patterns.benchmark_tier_cvt_filename(config)
    script:
        "../src/legendsimflow/scripts/tier/cvt.py"


rule plot_tier_cvt_observables:
    """Produces plots of observables from the tier `cvt`.

    Uses wildcard `simid`.
    """
    message:
        "Producing control plots for job cvt.{wildcards.simid}"
    input:
        rules.build_tier_cvt.output,
    output:
        patterns.plot_tier_cvt_observables_filename(config),
    script:
        "../src/legendsimflow/scripts/plots/tier_cvt_observables.py"
