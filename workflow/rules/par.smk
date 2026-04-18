# Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Rules to compute the simulation parameters (`par` step)."""

from legendsimflow import aggregate, hpge_pars, patterns


rule gen_all_tier_par:
    """Produce all `par` step outputs."""
    input:
        aggregate.gen_list_of_all_par_outputs(config),
        lambda wc: aggregate.gen_list_of_all_plots_outputs(
            config, tier="par", cache=smk_load_hpge_cache()
        ),


rule make_simstat_partition_file:
    """Create the simulation event statistics partitioning file.

    This rule maps chunks of event indices to partitions associated to the data
    taking runs specified in the "runlist" (from e.g. `config.runlist`) and
    stores them on disk as YAML files. The format is the following:

    ```yaml
    job_000:
      l200-p03-r001-phy: [0, 300]
      l200-p03-r002-phy: [301, 456]
    job_001:
      l200-p03-r002-phy: [0, 200]
      l200-p03-r003-phy: [201, 156]
    job_002:
      l200-p03-r003-phy: [0, 50]
    ```

    The events simulated in job `0` (456) are split between `r001` and `r002`.
    The partition corresponding to `r002` is however incomplete, and 200 events
    are taken from the simulation job `1`.

    The fraction of total simulated events (summed over all simulation jobs)
    that belong to a partition is determined by weighting with the fraction of
    livetime that belongs to that run.

    Uses wildcard `simid`.
    """
    input:
        stp_files=lambda wildcards: aggregate.gen_list_of_simid_outputs(
            config, tier="stp", simid=wildcards.simid
        ),
    params:
        # NOTE: these are not strictly needed here, but in this way Snakemake
        # can track these dependencies. these variables *do get used* in the
        # script.
        runinfo=config.metadata.datasets.runinfo,
        runlist=lambda wc: aggregate.get_runlist(config, wc.simid),
    output:
        patterns.simstat_part_filename(config),
    log:
        patterns.log_simstat_part_filename(config),
    script:
        "../src/legendsimflow/scripts/make_simstat_partition_file.py"


def smk_hpge_drift_time_map_inputs(wildcards):
    """Prepare inputs for the HPGe drift time map rule."""
    meta = config.metadata.hardware.detectors.germanium.diodes[wildcards.hpge_detector]
    ids = {"bege": "B", "coax": "C", "ppc": "P", "icpc": "V"}
    crystal_name = (
        ids[meta.type] + format(meta.production.order, "02d") + meta.production.crystal
    )

    _m = Path(config.paths.metadata)

    diode = _m / f"hardware/detectors/germanium/diodes/{wildcards.hpge_detector}.yaml"
    crystal = _m / f"hardware/detectors/germanium/crystals/{crystal_name}.yaml"

    return {
        "detdb_file": diode,
        "crydb_file": crystal,
        "dtmap_settings": _m
        / f"simprod/config/pars/{config.experiment}/geds/dtmap/settings.yaml",
        "_dummy": rules._init_julia_env.output,
    }


rule build_hpge_drift_time_map:
    """Produce an HPGe drift time map.

    Run a Julia script based on a pulse shape simulation performed with the
    `SolidStateDetectors.jl` package, using crystal geometry information from
    `legend-metadata`.

    Uses wildcards `hpge_detector` and `hpge_voltage`.
    """
    message:
        "Generating drift time map for HPGe detector {wildcards.hpge_detector} at {wildcards.hpge_voltage}V"
    input:
        unpack(smk_hpge_drift_time_map_inputs),
    params:
        metadata_path=config.paths.metadata,
    output:
        patterns.output_dtmap_filename(config),
    log:
        patterns.log_dtmap_filename(config),
    benchmark:
        patterns.benchmark_dtmap_filename(config)
    # NOTE: not using the `script` directive here since Snakemake has no nice
    # way to handle package dependencies nor Project.toml
    shell:
        "julia --project=workflow/src/LegendSimflow.jl --threads 1"
        "  workflow/src/legendsimflow/scripts/make_hpge_drift_time_maps.jl"
        "    --detector {wildcards.hpge_detector}"
        f"   --metadata {config.paths.metadata}"
        "    --dtmap-settings {input.dtmap_settings}"
        "    --opv {wildcards.hpge_voltage}"
        "    --output-file {output} &> {log}"


rule merge_hpge_drift_time_maps:
    """Merge HPGe drift time maps in a single file.

    Copy the top-level LH5 objects from each individual detector drift time map
    file into a single merged file using `h5copy`.

    Uses wildcard `runid`.
    """
    message:
        "Merging HPGe drift time map files for {wildcards.runid}"
    input:
        lambda wc: aggregate.gen_list_of_dtmaps(
            config, wc.runid, cache=smk_load_hpge_cache()
        ),
    output:
        patterns.output_dtmap_merged_filename(config),
    shell:
        r"""
        shopt -s nullglob
        out={output}

        # expand input files
        set -- {input}

        # if no matches, create an empty hdf5 file
        if [ "$#" -eq 0 ]; then
          python -c "import h5py; h5py.File('$out', 'w')"
          exit 0
        fi

        # seed with the first file
        cp "$1" "$out"
        shift

        # merge top-level objects from the rest
        for f in "$@"; do
          h5ls "$f" | awk '{{print $1}}' | while read -r o; do
            h5copy -i "$f" -o "$out" -s "/$o" -d "/$o"
          done
        done
        """


rule plot_hpge_drift_time_maps:
    """Produce a validation plot of an HPGe drift time map.

    Generates diagnostic plots of the computed drift time map for a single
    detector at the specified operational voltage.

    Uses wildcards `hpge_detector` and `hpge_voltage`.
    """
    message:
        "Plotting drift time map for HPGe {wildcards.hpge_detector} at {wildcards.hpge_voltage}V"
    input:
        rules.build_hpge_drift_time_map.output,
    output:
        patterns.plot_dtmap_filename(config),
    script:
        "../src/legendsimflow/scripts/plots/hpge_drift_time_maps.py"


rule extract_hpge_current_pulse_model:
    """Extract the HPGe current signal model for all detectors in a run.

    For each HPGe detector valid for modeling in the given run, either reads
    parameters from ``simprod/config/pars/geds/currmod/`` metadata (fast path,
    no l200data required) or fits current pulse waveforms from LEGEND-200 data.
    Results are written to a single YAML file keyed by detector name, and all
    fit validation plots are collected into one PDF.

    Uses wildcard ``runid``.

    Warning:
        when the script falls back to LEGEND-200 data, it discovers the
        underlying ``l200data`` input files dynamically at runtime. These files
        are therefore not tracked by Snakemake as rule inputs, so changes to
        them will not automatically trigger a rerun unless this rule is forced
        or those inputs are declared explicitly.
    """
    message:
        "Extracting current pulse model for all detectors in {wildcards.runid}"
    input:
        modelable_hpges=rules.cache_modelable_hpges.output[0],
    output:
        pars_file=patterns.output_currmod_merged_filename(config),
        plot_file=patterns.plot_currmod_filename(config),
    log:
        patterns.log_currmod_filename(config),
    script:
        "../src/legendsimflow/scripts/extract_hpge_current_pulse_model.py"


rule extract_hpge_observables_models:
    """Extract and store on disk models of the HPGe observables for a run.

    Stores YAML files with a mapping between HPGe detectors and respective
    information to reconstruct:

    - the energy resolution as a function of energy;
    - the A/E resolution as a function of energy;

    as determined during energy calibration. This is done in a separate rule
    because the data production parameter database is large and we don't want
    to use a lot of memory in the `build_tier_hit` rule.

    **Design:** this rule is a *collection* step, not a *validation* step. It
    gathers what it can from l200data and ``simprod/config/pars/geds/eresmod/``;
    the output may be incomplete. Completeness is validated downstream in
    ``build_tier_hit``.

    Uses wildcard `runid`.
    """
    message:
        "Extracting HPGe observables models for {wildcards.runid}"
    # NOTE: we don't list the file dependencies here because they are
    # dynamically generated, and that would slow down the DAG generation
    output:
        eresmod_file=patterns.output_eresmod_filename(config),
        aoeresmod_file=patterns.output_aoeresmod_filename(config),
        psdcuts_file=patterns.output_psdcuts_filename(config),
    log:
        patterns.log_eresmod_filename(config),
    script:
        "../src/legendsimflow/scripts/pars/extract_hpge_observables_models.py"
