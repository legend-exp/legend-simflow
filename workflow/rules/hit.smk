rule gen_all_tier_hit:
    """Aggregate and produce all the hit tier files."""
    input:
        aggregate.gen_list_of_all_simid_outputs(config, tier="hit"),


rule build_tier_hit:
    """Produces a hit tier file starting from a single `stp` tier file."""
    message:
        "Producing output file for job hit.{wildcards.simid}.{wildcards.jobid}"
    input:
        geom=patterns.geom_gdml_filename(config, tier="stp"),
        stp_file=patterns.output_simjob_filename(config, tier="stp"),
        optmap_lar=config.paths.optical_maps.lar,
        hpge_dtmaps=aggregate.gen_list_of_merged_dtmaps(config),
    output:
        patterns.output_simjob_filename(config, tier="hit"),
    log:
        patterns.log_filename(config, proctime, tier="hit"),
    benchmark:
        patterns.benchmark_filename(config, tier="hit")
    script:
        "../src/legendsimflow/scripts/tier/hit.py"


def smk_hpge_drift_time_map_inputs(wildcards):
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
    opv = _m / f"simprod/config/pars/opv/{runid_no_dt}-T%-all-opvs.yaml"

    return {
        "detdb_file": diode,
        "crydb_file": crystal,
        "opv_file": opv,
    }


rule build_hpge_drift_time_map:
    """Produce an HPGe drift time map.

    Uses wildcards `hpge_detector` and `runid`.
    """
    message:
        "Generating drift time map for HPGe detector {wildcards.hpge_detector} in run {wildcards.runid}"
    input:
        unpack(smk_hpge_drift_time_map_inputs),
    output:
        temp(patterns.output_dtmap_filename(config)),
    log:
        patterns.log_dtmap_filename(config, proctime),
    threads: 4
    params:
        metadata_path=config.paths.metadata,
    shell:
        "julia --project=workflow/src/legendsimflow/scripts --threads {threads}"
        "  workflow/src/legendsimflow/scripts/make_hpge_drift_time_maps.jl"
        "    --detector {wildcards.hpge_detector}"
        f"   --metadata {config.paths.metadata}"
        "    --opv-file {input.opv_file}"
        "    --output-file {output} &> {log}"


rule merge_hpge_drift_time_maps:
    """Merge HPGe drift time maps in a single file.

    Uses wildcard `runid`.
    """
    message:
        "Merging HPGe drift time map files for {wildcards.runid}"
    input:
        lambda wc: aggregate.gen_list_of_dtmaps(config, wc.runid),
    output:
        patterns.output_dtmap_merged_filename(config),
    params:
        input_regex=patterns.output_dtmap_filename(config, hpge_detector="*"),
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
          h5ls -1 "$f" | while read -r o; do
            h5copy -i "$f" -o "$out" -s "/$o" -d "/$o"
          done
        done
        """
