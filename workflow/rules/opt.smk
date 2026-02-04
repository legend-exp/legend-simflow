from dbetto.utils import load_dict

from legendsimflow import patterns, aggregate
from legendsimflow import metadata as mutils


rule gen_all_tier_opt:
    """Aggregate and produce all the opt tier files."""
    input:
        aggregate.gen_list_of_all_simid_outputs(config, tier="opt"),
        aggregate.gen_list_of_all_plots_outputs(config, tier="opt"),


# NOTE: we don't rely on rules from other tiers here (e.g.
# rules.build_tiers_stp.output) because we want to support making only the opt
# tier via the config.make_tiers option
rule build_tier_opt:
    """Produces a `opt` tier file starting from a single `stp` tier file.

    This rule implements the post-processing of the `stp` tier liquid argon
    energy depositions in chunks, in the following steps:

    - scintillation photons are generated corresponding to simulated energy
      depositions;
    - detected photoelectrons are sampled according to the input optical map;
    - a new time-coincidence map (TCM) across the processed SiPMs is created
      and stored in the output file.

    This rule can sample photoelectrons in each SiPM individually or for all
    SiPMs at the same time, see relevant `param` flag.

    The `stp` data format is preserved: SiPM tables are stored separately in
    the output file below `/hit/{sipm_name}`.

    Uses wildcards `simid` and `jobid`.
    """
    message:
        "Producing output file for job opt.{wildcards.simid}.{wildcards.jobid}"
    input:
        geom=patterns.geom_gdml_filename(config, tier="stp"),
        stp_file=patterns.output_simjob_filename(config, tier="stp"),
        optmap_lar=config.paths.optical_maps.lar,
    params:
        optmap_per_sipm=False,
        scintillator_volume_name="liquid_argon",
    output:
        patterns.output_simjob_filename(config, tier="opt"),
    log:
        patterns.log_filename(config, SIMFLOW_CONTEXT.proctime, tier="opt"),
    benchmark:
        patterns.benchmark_filename(config, tier="opt")
    script:
        "../src/legendsimflow/scripts/tier/opt.py"


rule plot_tier_opt_observables:
    """Produces plots of observables from the tier `opt`.

    Uses wildcard `simid`.
    """
    message:
        "Producing control plots for job opt.{wildcards.simid}"
    input:
        lambda wc: aggregate.gen_list_of_simid_outputs(
            config, tier="opt", simid=wc.simid
        ),
    output:
        patterns.plot_tier_opt_observables_filename(config),
    script:
        "../src/legendsimflow/scripts/plots/tier_opt_observables.py"
