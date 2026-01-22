import re

from dbetto.utils import load_dict

from legendsimflow import patterns, aggregate, hpge_pars
from legendsimflow import metadata as mutils
from legendsimflow import nersc


rule gen_all_tier_hit:
    """Aggregate and produce all the hit tier files."""
    input:
        aggregate.gen_list_of_all_simid_outputs(config, tier="hit"),
        aggregate.gen_list_of_all_plots_outputs(
            config, tier="hit", cache=SIMFLOW_CONTEXT.modelable_hpges
        ),


rule make_simstat_partition_file:
    """Creates the file containing the simulation event statistics partitioning file.

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
        config.paths.genmeta / "simstat" / "partitions_{simid}.yaml",
    log:
        patterns.log_simstat_part_filename(config, SIMFLOW_CONTEXT.proctime),
    script:
        "../src/legendsimflow/scripts/make_simstat_partition_file.py"


# NOTE: we don't rely on rules from other tiers here (e.g.
# rules.build_tiers_stp.output) because we want to support making only the hit
# tier via the config.make_tiers option
rule build_tier_hit:
    """Produce a `hit` tier file starting from a single `stp` tier file.

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
        # NOTE: technically this rule only depends on one block in the
        # partitioning file, but in practice the full file will always change
        simstat_part_file=rules.make_simstat_partition_file.output[0],
    output:
        patterns.output_simjob_filename(config, tier="hit"),
    log:
        patterns.log_filename(config, SIMFLOW_CONTEXT.proctime, tier="hit"),
    benchmark:
        patterns.benchmark_filename(config, tier="hit")
    script:
        "../src/legendsimflow/scripts/tier/hit.py"


def smk_hpge_drift_time_map_inputs(wildcards):
    """Prepare inputs for the HPGe drift time map rule."""
    meta = config.metadata.hardware.detectors.germanium.diodes[wildcards.hpge_detector]
    ids = {"bege": "B", "coax": "C", "ppc": "P", "icpc": "V"}
    crystal_name = (
        ids[meta.type] + format(meta.production.order, "02d") + meta.production.crystal
    )

    # remove the datatype at the end of the runid string, it's not needed to
    # locate the operational voltage file
    runid_no_dt = "-".join(wildcards.runid.split("-")[:-1])

    _m = Path(config.paths.metadata)

    diode = (
        _m / f"hardware/detectors/germanium/diodes/{wildcards.hpge_detector}.yaml",
    )
    crystal = _m / f"hardware/detectors/germanium/crystals/{crystal_name}.yaml"

    return {
        "detdb_file": diode,
        "crydb_file": crystal,
        "_dummy": rules._init_julia_env.output,
    }


rule build_hpge_drift_time_map:
    """Produce an HPGe drift time map.

    Uses wildcards `hpge_detector` and `runid`.
    """
    message:
        "Generating drift time map for HPGe detector {wildcards.hpge_detector} in run {wildcards.runid}"
    input:
        unpack(smk_hpge_drift_time_map_inputs),
    params:
        metadata_path=config.paths.metadata,
        operational_voltage=lambda wc: mutils.simpars(
            config.metadata, "geds.opv", wc.runid
        )[wc.hpge_detector].operational_voltage_in_V,
    output:
        temp(patterns.output_dtmap_filename(config)),
    log:
        patterns.log_dtmap_filename(config, SIMFLOW_CONTEXT.proctime),
    benchmark:
        patterns.benchmark_dtmap_filename(config)
    # NOTE: not using the `script` directive here since Snakemake has no nice
    # way to handle package dependencies nor Project.toml
    shell:
        "julia --project=workflow/src/legendsimflow/scripts --threads 1"
        "  workflow/src/legendsimflow/scripts/make_hpge_drift_time_maps.jl"
        "    --detector {wildcards.hpge_detector}"
        f"   --metadata {config.paths.metadata}"
        "    --opv {params.operational_voltage}"
        "    --output-file {output} &> {log}"


rule plot_hpge_drift_time_maps:
    """Produce a validation plot of an HPGe drift time map.

    Uses wildcards `hpge_detector` and `runid`.
    """
    message:
        "Plotting HPGe drift time map for {wildcards.hpge_detector} in {wildcards.runid}"
    input:
        rules.build_hpge_drift_time_map.output,
    output:
        patterns.plot_dtmap_filename(config),
    script:
        "../src/legendsimflow/scripts/plots/hpge_drift_time_maps.py"


rule merge_hpge_drift_time_maps:
    """Merge HPGe drift time maps in a single file.

    Uses wildcard `runid`.
    """
    message:
        "Merging HPGe drift time map files for {wildcards.runid}"
    input:
        lambda wc: aggregate.gen_list_of_dtmaps(
            config, wc.runid, cache=SIMFLOW_CONTEXT.modelable_hpges
        ),
    params:
        input_regex=lambda wc: patterns.output_dtmap_filename(
            config, runid=wc.runid, hpge_detector="*"
        ),
    output:
        patterns.output_dtmap_merged_filename(config),
    shell:
        r"""
        shopt -s nullglob
        out={output}

        # expand glob into $1 $2 ...
        set -- {params.input_regex}

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

    with evt_idx_file.open() as f:
        f.write(wf_idx)

    return {
        "raw_file": raw_file,
        "raw_wf_idx_file": evt_idx_file,
        "dsp_cfg_file": dsp_cfg_file,
    }


rule extract_current_pulse_model:
    """Extract the HPGe current signal model.

    Perform a fit of the current waveform and stores the best-fit model
    parameters in a YAML file.

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
        patterns.log_currmod_filename(config, SIMFLOW_CONTEXT.proctime),
    script:
        "../src/legendsimflow/scripts/extract_hpge_current_pulse_model.py"


rule merge_current_pulse_model_pars:
    """Merge the HPGe current signal model parameters in a single file per `runid`.

    Uses wildcard `runid`.
    """
    message:
        "Merging current model parameters in {wildcards.runid}"
    input:
        lambda wc: aggregate.gen_list_of_currmods(
            config, wc.runid, cache=SIMFLOW_CONTEXT.modelable_hpges
        ),
    output:
        patterns.output_currmod_merged_filename(config),
    run:
        import dbetto

        # NOTE: this is guaranteed to be sorted as in the input file list
        hpges = SIMFLOW_CONTEXT.modelable_hpges[wildcards.runid]

        out_dict = {}
        for i, f in enumerate(input):
            out_dict[hpges[i]] = dbetto.utils.load_dict(f)

        dbetto.utils.write_dict(out_dict, output[0])
