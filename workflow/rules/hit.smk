from legendsimflow import aggregate, patterns
from legendsimflow.metadata import get_tier_settings

# realistic-PSL-based PSD simulation keeps a detector's PSL waveforms in
# memory, pushing peak RSS above the usual 2 GB budget. left empty otherwise,
# so default-resources applies
_hit_resources = (
    {"mem_mb": 2300}
    if get_tier_settings(config, "hit").get("simulate_psd_with_psl", False)
    else {}
)

# gates the whole pre-correction -> fit -> merge -> final-hit-input chain (see
# smk_build_tier_hit_extra_input below). Off by default: the final hit is built
# once, with per-detector corrections falling back to the aoemeanmod_default in
# the hit tier settings.
_two_pass_aoe_corr = get_tier_settings(config, "hit").get(
    "two_pass_aoe_correction", False
)


def smk_hit_inputs(wc):
    """Inputs shared by both hit-tier builds (pass-1 pre-correction and pass-2 final).

    Factored out so the pre-correction rule can reuse them verbatim via
    ``use rule ... with:``. The two builds differ only in output path and, for
    the final build, the optional ``aoemean_file`` input.
    """
    hit_settings = get_tier_settings(config, "hit")
    return {
        "geom": patterns.geom_gdml_filename(config, tier="stp", simid=wc.simid),
        "stp_file": patterns.output_simjob_filename(
            config, tier="stp", simid=wc.simid, jobid=wc.jobid
        ),
        "hpge_dtmaps": aggregate.gen_list_of_merged_dtmaps(
            config, wc.simid, has_psd=hit_settings.get("simulate_psd", True)
        ),
        "hpge_currmods": aggregate.gen_list_of_merged_currmods(
            config, wc.simid, has_psd=hit_settings.get("simulate_psd", True)
        ),
        "hpge_eresmods": aggregate.gen_list_of_eresmods(config, wc.simid),
        "hpge_aoeresmods": aggregate.gen_list_of_aoeresmods(config, wc.simid),
        "hpges_realistic_psls": aggregate.gen_list_of_merged_realistic_psls(
            config,
            wc.simid,
            simulate_psd_with_psl=hit_settings.get("simulate_psd_with_psl", False),
        ),
        "hpge_psdcuts": aggregate.gen_list_of_psdcuts(config, wc.simid),
        # NOTE: technically this rule only depends on one block in the
        # partitioning file, but in practice the full file will always change
        "simstat_part_file": patterns.simstat_part_filename(config, simid=wc.simid),
        "usability": rules.cache_detector_usabilities.output.usability,
        "psd_usability": rules.cache_detector_usabilities.output.psd_usability,
        "crystal_metadata_usability": rules.cache_detector_usabilities.output.crystal_metadata_usability,
    }


def smk_build_tier_hit_extra_input(wc):
    """The final hit build's only extra input: the merged A/E correction, toggle-gated.

    Returning ``{}`` (not ``{"aoemean_file": []}``) when off is load-bearing: an
    absent named input is mapped to ``None`` by ``snakemake_argparse_bridge``,
    whereas an empty-list value would stringify to the literal ``"[]"``.
    """
    if _two_pass_aoe_corr:
        return {"aoemean_file": patterns.output_aoemeanmod_merged_filename(config)}
    return {}


rule gen_all_tier_hit:
    """Aggregate and produce all the hit tier files."""
    input:
        aggregate.gen_list_of_all_simid_outputs(config, tier="hit"),
        lambda wc: aggregate.gen_list_of_all_plots_outputs(config, tier="hit"),
        # SSD-modeling provenance detinfo, aggregated from the drift-time-map
        # sidecars; gated by the same `simulate_psd` setting that builds them
        lambda wc: (
            [patterns.detinfo_filename(config, "hpge_ssd_modeling")]
            if get_tier_settings(config, "hit").get("simulate_psd", True)
            else []
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

    When the `two_pass_aoe_correction` hit-tier setting is enabled, this final
    build additionally consumes the merged A/E energy-dependence correction
    (`merge_aoe_energy_corrections`), computed from a pre-correction pass (see
    `build_tier_hit_precorr`).

    Uses wildcards `simid` and `jobid`.
    """
    message:
        "Producing output file for job hit.{wildcards.simid}.{wildcards.jobid}"
    input:
        unpack(smk_hit_inputs),
        unpack(smk_build_tier_hit_extra_input),
    output:
        patterns.output_simjob_filename(config, tier="hit"),
    log:
        patterns.log_filename(config, tier="hit"),
    benchmark:
        patterns.benchmark_filename(config, tier="hit")
    resources:
        **_hit_resources,
    script:
        "../src/legendsimflow/scripts/tier/hit.py"


# pass 1 of the two-pass A/E energy correction: identical post-processing to
# build_tier_hit but with no correction applied (input overridden to drop the
# aoemean_file, which also breaks the would-be dependency cycle) and a temporary
# output. Only enters the DAG when two_pass_aoe_correction is enabled, i.e. when
# extract_hpge_aoemean_energy_dependence requests these files.
use rule build_tier_hit as build_tier_hit_precorr with:
    message:
        "Producing pre-correction hit output for job hit.{wildcards.simid}.{wildcards.jobid}"
    input:
        unpack(smk_hit_inputs),
    output:
        temp(patterns.output_simjob_precorr_hit_filename(config)),
    log:
        patterns.log_filename(config, tier="hit_precorr"),
    benchmark:
        patterns.benchmark_filename(config, tier="hit_precorr")


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
