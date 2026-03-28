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


def smk_extract_current_pulse_model_inputs(wildcards):
    """Prepare inputs for the HPGe current model extraction rule."""
    raw_file, wf_idx, dsp_cfg_file = hpge_pars.find_current_pulse_model_inputs(
        config,
        wildcards.runid,
        wildcards.hpge_detector,
        hit_tier_name="hit",
        use_hpge_name="true",
    )

    evt_idx_file = patterns.input_currmod_evt_idx_file(
        config, runid=wildcards.runid, hpge_detector=wildcards.hpge_detector
    )

    with evt_idx_file.open("w") as f:
        f.write(wf_idx)

    return {
        "raw_file": raw_file,
        "raw_wf_idx_file": evt_idx_file,
        "dsp_cfg_file": dsp_cfg_file,
    }


rule extract_current_pulse_model:
    """Extract the HPGe current signal model.

    Perform a fit of current signals recorded in LEGEND-200 and stores the
    best-fit model parameters in a YAML file.

    :::{warning}
    This rule does not have the relevant LEGEND-200 data files as input, since
    they are dynamically discovered and this would therefore slow down the DAG
    generation. Therefore, remember to force-rerun if the input data is
    updated!
    :::

    Uses wildcards `runid` and `hpge_detector`.
    """
    message:
        "Extracting current model for detector {wildcards.hpge_detector} in {wildcards.runid}"
    # NOTE: we don't list the file dependencies here because they are
    # dynamically generated, and that would slow down the DAG generation
    # input:
    #     unpack(smk_extract_current_pulse_model_inputs),
    output:
        pars_file=temp(patterns.output_currmod_filename(config)),
        plot_file=patterns.plot_currmod_filename(config),
    log:
        patterns.log_currmod_filename(config),
    script:
        "../src/legendsimflow/scripts/extract_hpge_current_pulse_model.py"


rule merge_current_pulse_model_pars:
    """Merge the HPGe current signal model parameters in a single file per `runid`.

    Collect the individual best-fit parameter files (one per detector) and
    write them into a single YAML file keyed by detector name.

    Uses wildcard `runid`.
    """
    message:
        "Merging current model parameters in {wildcards.runid}"
    input:
        lambda wc: aggregate.gen_list_of_currmods(
            config, wc.runid, cache=smk_load_hpge_cache()
        ),
    params:
        # materialize the HPGe list here so the run: block can map file
        # indices back to detector names without calling back into the cache
        hpges=lambda wc: list(smk_load_hpge_cache()[wc.runid]),
    output:
        patterns.output_currmod_merged_filename(config),
    run:
        import dbetto

        # NOTE: this is guaranteed to be sorted as in the input file list
        hpges = params.hpges

        out_dict = {}
        for i, f in enumerate(input):
            out_dict[hpges[i]] = dbetto.utils.load_dict(f)

        dbetto.utils.write_dict(out_dict, output[0])


rule extract_hpge_observables_models:
    """Extract from the LEGEND-200 data and store on disk models of the HPGe observables.

    Stores YAML files with a mapping between HPGe detectors and respective
    information to reconstruct:

    - the energy resolution as a function of energy;
    - the A/E resolution as a function of energy;

    as determined during energy calibration. This is done in a separate rule
    because the data production parameter database is large and we don't want
    to use a lot of memory in the `build_tier_hit` rule.

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
    run:
        import dbetto
        from legendsimflow import hpge_pars, utils

        l200data = config.paths.l200data

        hit_tier_name = utils.get_hit_tier_name(l200data)
        pars_db = utils.init_generated_pars_db(l200data, tier=hit_tier_name, lazy=True)

        eres_pars_dict = hpge_pars.lookup_energy_res_metadata(
            l200data,
            config.metadata,
            wildcards.runid,
            hit_tier_name=hit_tier_name,
            pars_db=pars_db,
        )

        out_dict = dbetto.AttrsDict({})
        fields = ["expression", "parameters", "uncertainties"]
        for hpge, meta in eres_pars_dict.items():
            out_dict[hpge] = {f: meta[f] for f in fields}

        dbetto.utils.write_dict(out_dict.to_dict(), output.eresmod_file)

        aoeres_pars_dict = hpge_pars.lookup_aoe_res_metadata(
            l200data,
            config.metadata,
            wildcards.runid,
            hit_tier_name=hit_tier_name,
            pars_db=pars_db,
        )

        out_dict = dbetto.AttrsDict({})
        fields = ["expression", "pars", "errs"]
        for hpge, meta in aoeres_pars_dict.items():
            out_dict[hpge] = {f: meta[f] for f in fields}

        dbetto.utils.write_dict(out_dict.to_dict(), output.aoeresmod_file)

        aoecuts_pars_dict = hpge_pars.lookup_psd_cut_values(
            l200data,
            config.metadata,
            wildcards.runid,
            hit_tier_name=hit_tier_name,
            pars_db=pars_db,
        )

        dbetto.utils.write_dict(aoecuts_pars_dict, output.psdcuts_file)
