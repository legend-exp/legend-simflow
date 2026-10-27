from collections.abc import Mapping
from functools import partial

from dbetto.utils import load_dict

from legendsimflow import patterns, aggregate
from legendsimflow import metadata as mutils
from legendsimflow.exceptions import SimflowConfigError
from legendsimflow.metadata import (
    deferred_tier_setting,
    get_opt_scintillators,
    get_tier_settings,
)

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


def smk_optmaps(wildcards):
    """The optical maps of the light sources, ``optmap_<source> -> file``.

    Not needed when the light comes from tracked photons, nor for the sources
    with a light model.
    """
    if get_tier_settings(config, "opt").light_source == "tracked_photons":
        return {}

    sources = {
        name: source
        for name, source in get_opt_scintillators(config).items()
        if source.model is None
    }
    if "lar_patch" in config.paths.optical_maps and "lar" not in sources:
        msg = (
            "paths.optical_maps.lar_patch patches the map of the lar light source, "
            f"but no light source with an optical map has this name: {sorted(sources)}"
        )
        raise SimflowConfigError(msg)

    optmaps = {}
    for name in sources:
        if name not in config.paths.optical_maps:
            msg = f"paths.optical_maps.{name} is needed by the {name} light source"
            raise SimflowConfigError(msg)
        if name == "lar":
            optmaps["optmap_lar"] = _optmap_lar(config, wildcards.simid)
        else:
            optmaps[f"optmap_{name}"] = on_scratch_smk(config.paths.optical_maps[name])
    return optmaps


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
        # read only when the rule is used: there might be no lar map
        base=lambda wc: on_scratch_smk(config.paths.optical_maps.lar),
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

    This rule implements the post-processing of the `stp` tier energy
    depositions in the scintillators (liquid argon, PEN, ...) in chunks, in the
    following steps:

    - each chunk is partitioned according to the livetime span of each run
      (see the `make_simstat_partition_file` rule). For each partition:
    - the detector usability is retrieved from `legend-metadata` and stored in
      the output;
    - with `light_source: optmap` (default), scintillation photons are generated
      corresponding to simulated energy depositions and detected
      photoelectrons are sampled according to the optical map of each light
      source (see `scintillators` in the tier settings). A light model (a
      Python function giving the photon detection probability) can replace
      the map of a light source, which then needs no map file. The
      photoelectrons of all sources in the same `stp` event window (TCM row)
      are pooled, and each window gives one output row;
    - with `light_source: tracked_photons`, the photons recorded in the `stp`
      SiPM tables are used instead and the optical map is not an input. As in
      the `hit` tier, each SiPM row gives one output row, and every photon is
      a photoelectron: the SiPM efficiency is applied while tracking. The
      photoelectrons of each SiPM row are capped as with the map. SiPMs that
      no photon reached get an empty table;
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
        unpack(smk_optmaps),
        geom=patterns.geom_gdml_filename(config, tier="stp"),
        # read by the light models that depend on the geometry layout
        geom_config=patterns.geom_config_filename(config, tier="stp"),
        stp_file=patterns.output_simjob_filename(config, tier="stp"),
        # NOTE: technically this rule only depends on one block in the
        # partitioning file, but in practice the full file will always change
        simstat_part_file=patterns.simstat_part_filename(config),
        usability=rules.cache_detector_usabilities.output.usability,
    params:
        optmap_per_sipm=_tier_setting("opt", "optmap_per_sipm"),
        # holds the optmap_scaling_factor of each light source
        scintillators=lambda wc: get_opt_scintillators(config),
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
