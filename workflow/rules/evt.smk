from legendsimflow import patterns, aggregate


rule gen_all_tier_evt:
    """Aggregate and produce all the evt tier files."""
    input:
        aggregate.gen_list_of_all_simid_outputs(config, tier="evt"),
        # aggregate.gen_list_of_all_tier_evt_outputs(config),


rule build_tier_evt:
    """Produces an evt tier file."""
    message:
        "Producing output file for job evt.{wildcards.simid}.{wildcards.jobid}"
    input:
        opt_file=patterns.output_simjob_filename(config, tier="opt"),
        hit_file=patterns.output_simjob_filename(config, tier="hit"),
    params:
        buffer_len="500*MB",
    output:
        patterns.output_simjob_filename(config, tier="evt"),
    log:
        patterns.log_filename(config, SIMFLOW_CONTEXT.proctime, tier="evt"),
    benchmark:
        patterns.benchmark_filename(config, tier="evt")
    script:
        "../src/legendsimflow/scripts/tier/evt.py"


# rule build_concat_tier_evt:
#     """Produces a concatenated evt tier file."""
#     message:
#         "Producing output file for job evt.{wildcards.simid}"
#     input:
#         evt_files=lambda wildcards: aggregate.gen_list_of_simid_outputs(
#             config, tier="evt", simid=wildcards.simid
#         ),
#     params:
#         buffer_len="500*MB",
#     output:
#         patterns.output_evt_filename(config),
#     log:
#         patterns.log_evtfile_path(config, proctime),
#     benchmark:
#         patterns.benchmark_evtfile_path(config)
#     script:
#         "boh"
