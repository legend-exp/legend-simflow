from collections.abc import Mapping
from functools import partial

from dbetto.utils import load_dict

from legendsimflow import patterns, aggregate
from legendsimflow import metadata as mutils
from legendsimflow.metadata import deferred_tier_setting

_tier_setting = partial(deferred_tier_setting, config)


def _optmap_patch(config, simid):
    # a single map for every simid, or a mapping <simid> -> <patch map>; no patch is []
    patch = config.paths.optical_maps.get("lar_patch", [])

    if not isinstance(patch, Mapping):
        return patch

    return patch.get(simid, [])


def _optmap_lar(config, simid):
    # the patched map is already on scratch, if enabled
    if _optmap_patch(config, simid):
        return patterns.patched_optmap_filename(config, simid=simid)

    return on_scratch_smk(config.paths.optical_maps.lar)


rule gen_all_tier_opt:
    """Aggregate and produce all the opt tier files."""
    input:
        aggregate.gen_list_of_all_simid_outputs(config, tier="opt"),
        aggregate.gen_list_of_all_plots_outputs(config, tier="opt"),


rule patch_optical_map:
    """Substitute a separately simulated region into the LAr optical map.

    The base map is simulated without hardware that is only present in some runs
    -- a calibration source and its absorber, say -- so its detection
    probabilities are wrong in the volume around it. This rule replaces that
    region with a map of a smaller volume simulated with the hardware in place.

    Only runs for simids that configure a patch in
    ``paths.optical_maps.lar_patch``; the opt tier of the others reads the base
    map directly.

    Uses wildcard `simid`.
    """
    message:
        "Patching the LAr optical map for {wildcards.simid}"
    input:
        base=on_scratch_smk(config.paths.optical_maps.lar),
        patch=lambda wc: on_scratch_smk(_optmap_patch(config, wc.simid)),
    output:
        temp(patterns.patched_optmap_filename(config)),
    log:
        patterns.patched_optmap_log_filename(config),
    shell:
        "reboost-optical -v patchmap {input.base} {input.patch} {output} &> {log}"


# NOTE: we don't rely on rules from other tiers here (e.g.
# rules.build_tiers_stp.output) because we want to support making only the opt
# tier via the config.make_steps option
rule build_tier_opt:
    """Produce a `opt` tier file starting from a single `stp` tier file.

    This rule implements the post-processing of the `stp` tier liquid argon
    energy depositions in chunks, in the following steps:

    - each chunk is partitioned according to the livetime span of each run
      (see the `make_simstat_partition_file` rule). For each partition:
    - the detector usability is retrieved from `legend-metadata` and stored in
      the output;
    - scintillation photons are generated corresponding to simulated energy
      depositions;
    - detected photoelectrons are sampled according to the input optical map;
    - a finite resolution is applied to each photoelectron amplitude (see
      script);
    - photoelectrons are clustered in time to simulate the effect of finite
      time resolution of the system;
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
        optmap_lar=lambda wc: _optmap_lar(config, wc.simid),
        # NOTE: technically this rule only depends on one block in the
        # partitioning file, but in practice the full file will always change
        simstat_part_file=patterns.simstat_part_filename(config),
        usability=rules.cache_detector_usabilities.output.usability,
    params:
        optmap_per_sipm=_tier_setting("opt", "optmap_per_sipm"),
        scintillator_volume_name=_tier_setting("opt", "scintillator_volume_name"),
        optmap_scaling_factor=_tier_setting("opt", "optmap_scaling_factor"),
        photoelectron_resolution_sigma=_tier_setting(
            "opt", "photoelectron_resolution_sigma"
        ),
        time_resolution_in_ns=_tier_setting("opt", "time_resolution_in_ns"),
        max_pes_per_hit_per_sipm=_tier_setting("opt", "max_pes_per_hit_per_sipm"),
        max_pes_per_hit_combined=_tier_setting("opt", "max_pes_per_hit_combined"),
    output:
        patterns.output_simjob_filename(config, tier="opt"),
    log:
        patterns.log_filename(config, tier="opt"),
    benchmark:
        patterns.benchmark_filename(config, tier="opt")
    script:
        "../src/legendsimflow/scripts/tier/opt.py"


rule plot_tier_opt_observables:
    """Produce validation plots of observable distributions from the `opt` tier.

    Generates diagnostic plots from all `opt` output files for the given
    `simid`.

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
