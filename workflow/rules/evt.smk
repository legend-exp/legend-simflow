from legendsimflow import patterns, aggregate
from legendsimflow.metadata import get_tier_settings

_evt_settings = get_tier_settings(config, "evt")
_skip_opt = _evt_settings.get("skip_opt", False)
_skip_hit = _evt_settings.get("skip_hit", False)

if _skip_opt and _skip_hit:
    raise WorkflowError("evt: skip_opt and skip_hit cannot both be True")


def _tier_setting(tier, key):
    return lambda wc: config.metadata.simprod.config.tier[tier][
        config.experiment
    ].settings[key]


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

    A top-level `detector_uids` struct mapping detector names to reboost UIDs
    (union of hit and opt tiers) is also written, to enable downstream
    per-group filtering in the `pdf` tier.

    Uses wildcards `simid` and `jobid`.
    """
    message:
        "Producing output file for job evt.{wildcards.simid}.{wildcards.jobid}"
    input:
        stp_file=patterns.output_simjob_filename(config, tier="stp"),
        opt_file=lambda wc: (
            []
            if _skip_opt
            else [
                patterns.output_simjob_filename(
                    config, tier="opt", simid=wc.simid, jobid=wc.jobid
                )
            ]
        ),
        hit_file=lambda wc: (
            []
            if _skip_hit
            else [
                patterns.output_simjob_filename(
                    config, tier="hit", simid=wc.simid, jobid=wc.jobid
                )
            ]
        ),
        simstat_part_file=patterns.simstat_part_filename(config),
        usability=rules.cache_detector_usabilities.output.usability,
    params:
        add_random_coincidences=_tier_setting("evt", "add_random_coincidences"),
        geds_energy_thr_kev=_tier_setting("evt", "geds_energy_thr_kev"),
        spms_energy_thr_pe=_tier_setting("evt", "spms_energy_thr_pe"),
        skip_opt=_skip_opt,
        skip_hit=_skip_hit,
        # with add_random_coincidences, this rule reads forced-trigger evt
        # files discovered dynamically from l200data and not listed as inputs;
        # track the path so it reruns on change (None otherwise to avoid
        # spurious reruns)
        _l200data=(
            config.paths.get("l200data", None)
            if _evt_settings.get("add_random_coincidences", False)
            else None
        ),
    output:
        patterns.output_simjob_filename(config, tier="evt"),
    log:
        patterns.log_filename(config, tier="evt"),
    benchmark:
        patterns.benchmark_filename(config, tier="evt")
    script:
        "../src/legendsimflow/scripts/tier/evt.py"
