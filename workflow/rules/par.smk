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

from pathlib import Path

from legendsimflow import aggregate, hpge_pars, patterns, utils
from legendsimflow.metadata import (
    electron_gun_macro_template,
    electron_gun_primaries,
    get_par_settings,
    get_tier_settings,
)


rule gen_all_tier_par:
    """Produce all `par` step outputs."""
    input:
        aggregate.gen_list_of_all_par_outputs(config),
        lambda wc: aggregate.gen_list_of_all_plots_outputs(config, tier="par"),


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


def smk_hpge_psd_simulation_inputs(wildcards):
    """Prepare inputs for the HPGe PSD simulation (SSD) rules."""
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
        "ssd_settings": _m
        / f"simprod/config/pars/{config.experiment}/geds/ssd/settings.yaml",
        "_dummy": rules._init_julia_env.output,
    }


rule build_hpge_drift_time_map:
    """Produce an HPGe drift-time map.

    Run a Julia script based on a pulse-shape simulation performed with the
    [`SolidStateDetectors.jl`](https://juliaphysics.github.io/SolidStateDetectors.jl/stable/)
    package, using crystal geometry information from `legend-metadata`. The drift time recorded at each grid point is the time at
    which the simulated hole charge-collection signal reaches its maximum current
    (the argmax of the charge-signal derivative); one map is produced per
    crystal-axis angle.

    The output is one LH5 file per `(detector, voltage)` pair with a single
    top-level group named after the detector, holding the gridded map in the
    `(r, z)`-field format read back by
    {func}`reboost.hpge.utils.get_hpge_rz_field` and consumed by
    {func}`reboost.hpge.psd.drift_time` when the `hit` tier is built:

    | Field                    | Type         | Units | Description                                                                                                       |
    | ------------------------ | ------------ | ----- | --------------------------------------------------------------------------------------------------------------- |
    | `r`                      | `Array`      | m     | Radial coordinates of the rectangular `(r, z)` grid.                                                             |
    | `z`                      | `Array`      | m     | Axial coordinates of the grid. The origin `(r, z) = (0, 0)` is the centre of the p+ contact.                     |
    | `drift_time_<angle>_deg` | `Array` (2D) | ns    | Drift time over the `(r, z)` grid at crystal-axis azimuth `<angle>`. Pixels outside the detector profile hold NaN. |

    Uses wildcards `hpge_detector` and `hpge_voltage`.
    """
    message:
        "Generating drift-time map for HPGe detector {wildcards.hpge_detector} at {wildcards.hpge_voltage}V"
    input:
        unpack(smk_hpge_psd_simulation_inputs),
    params:
        metadata_path=config.paths.metadata,
    output:
        dtmap_file=patterns.output_dtmap_filename(config),
        info_file=patterns.output_dtmap_info_filename(config),
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
        "    --ssd-settings {input.ssd_settings}"
        "    --opv {wildcards.hpge_voltage}"
        "    --output-file {output.dtmap_file}"
        "    --info-file {output.info_file} &> {log}"


rule merge_hpge_drift_time_maps:
    """Merge HPGe drift-time maps in a single file.

    Copy the top-level LH5 objects from each individual detector drift-time map
    file into a single merged file using `h5copy`.

    Uses wildcard `runid`.
    """
    message:
        "Merging HPGe drift-time map files for {wildcards.runid}"
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


rule aggregate_hpge_ssd_modeling_info:
    """Aggregate the HPGe SSD-modeling provenance into a `detinfo` file.

    Each `build_hpge_drift_time_map` job writes, next to its LH5 map, a small
    YAML sidecar holding the SSD simulation provenance scalars for that
    `(detector, voltage)` pair (impurity scaling factor, measured and simulated
    depletion voltages). This rule collects those sidecars and pivots them into
    a single `pars/detinfo/hpge_ssd_modeling.yaml` mapping `runid -> detector ->
    {impurity_scaling_factor, measured_depletion_voltage_in_V,
    simulated_depletion_voltage_raw_in_V, simulated_depletion_voltage_in_V}`,
    following the `detinfo` convention (see {ref}`par-detinfo`).

    No wildcards are used.
    """
    localrule: True
    message:
        "Aggregating HPGe SSD-modeling provenance into a detinfo file"
    input:
        lambda wc: aggregate.gen_list_of_dtmap_info_files(
            config, cache=smk_load_hpge_cache()
        ),
    output:
        patterns.detinfo_filename(config, "hpge_ssd_modeling"),
    run:
        import dbetto

        info = aggregate.collect_hpge_ssd_modeling_info(config, smk_load_hpge_cache())
        dbetto.utils.write_dict(info, output[0])


rule plot_hpge_drift_time_maps:
    """Produce a validation plot of an HPGe drift-time map.

    Generates diagnostic plots of the computed drift-time map for a single
    detector at the specified operational voltage.

    Uses wildcards `hpge_detector` and `hpge_voltage`.
    """
    message:
        "Plotting drift-time map for HPGe {wildcards.hpge_detector} at {wildcards.hpge_voltage}V"
    input:
        rules.build_hpge_drift_time_map.output.dtmap_file,
    output:
        patterns.plot_dtmap_filename(config),
    script:
        "../src/legendsimflow/scripts/plots/hpge_drift_time_maps.py"


rule build_hpge_pulse_shape_library:
    """Produce an (ideal) HPGe pulse-shape library.

    Run a Julia script based on a pulse-shape simulation performed with the
    [`SolidStateDetectors.jl`](https://juliaphysics.github.io/SolidStateDetectors.jl/stable/)
    package, using crystal geometry information from `legend-metadata`. It simulates the full charge-collection transient at every
    `(r, z, angle)` grid point. These are un-convolved (ideal) charge waveforms;
    the electronics transfer function is applied later, in
    `convolve_hpge_ideal_pulse_shape_lib`.

    The output is one LH5 file per `(detector, voltage)` pair with a single
    top-level group named after the detector:

    | Field                  | Type         | Units | Description                                                                  |
    | ---------------------- | ------------ | ----- | --------------------------------------------------------------------------- |
    | `r`                    | `Array`      | m     | Radial coordinates of the simulation grid.                                  |
    | `z`                    | `Array`      | m     | Axial coordinates of the simulation grid.                                   |
    | `dt`                   | `Scalar`     | ns    | Time step between waveform samples.                                         |
    | `waveform_<angle>_deg` | `Array` (3D) |       | Ideal charge waveforms at azimuth `<angle>`, shape `(n_z, n_r, n_samples)`. |

    The sample times of a waveform are `t_i = i * dt` (no `t0` offset is stored).

    Uses wildcards `hpge_detector` and `hpge_voltage`.
    """
    message:
        "Generating pulse-shape library for HPGe detector {wildcards.hpge_detector} at {wildcards.hpge_voltage}V"
    input:
        unpack(smk_hpge_psd_simulation_inputs),
    params:
        metadata_path=config.paths.metadata,
    output:
        patterns.output_ideal_psl_filename(config),
    log:
        patterns.log_ideal_psl_filename(config),
    benchmark:
        patterns.benchmark_ideal_psl_filename(config)
    # NOTE: not using the `script` directive here since Snakemake has no nice
    # way to handle package dependencies nor Project.toml
    shell:
        "julia --project=workflow/src/LegendSimflow.jl --threads 1"
        "  workflow/src/legendsimflow/scripts/make_hpge_ideal_pulse_shape_lib.jl"
        "    --detector {wildcards.hpge_detector}"
        f"   --metadata {config.paths.metadata}"
        "    --ssd-settings {input.ssd_settings}"
        "    --opv {wildcards.hpge_voltage}"
        "    --output-file {output} &> {log}"


rule convolve_hpge_ideal_pulse_shape_lib:
    """Convolve ideal HPGe pulse-shape libraries with electronics response.

    Produce one realistic pulse-shape library per detector and run by selecting
    the ideal library for the detector operational voltage in that run and
    convolving it with the electronics parameters (`sigma`, `tau`) from a merged
    electronics-model YAML file, then differentiating, moving-window-averaging,
    aligning, and normalising to obtain current waveforms.

    The output is one LH5 file per `(runid, detector)` pair with a top-level
    group named after the detector, in the pulse-shape-library format read back
    by {func}`reboost.hpge.utils.get_hpge_pulse_shape_library` when the `hit`
    tier is built:

    | Field                    | Type         | Units | Description                                                                                            |
    | ------------------------ | ------------ | ----- | ----------------------------------------------------------------------------------------------------- |
    | `r`                      | `Array`      | mm    | Radial coordinates of the simulation grid (converted from metres).                                    |
    | `z`                      | `Array`      | mm    | Axial coordinates of the simulation grid (converted from metres).                                     |
    | `dt`                     | `Scalar`     | ns    | Sampling time step of the realistic current waveforms.                                                |
    | `t0`                     | `Scalar`     | ns    | Global time offset: time of the alignment index relative to the current peak.                         |
    | `waveform_<angle>_deg`   | `Array` (3D) |       | Processed, normalised current waveforms at azimuth `<angle>`, shape `(n_r, n_z, n_samples)`; NaN outside the detector. |
    | `drift_time_<angle>_deg` | `Array` (2D) | ns    | Drift time per `(r, z)` pixel at azimuth `<angle>`, shape `(n_r, n_z)`; NaN outside the detector.       |

    The current-waveform sample times are `t_i = t0 + i * dt`.

    Uses wildcards `runid` and `hpge_detector`.
    """
    message:
        "Convolving ideal pulse-shape library for detector {wildcards.hpge_detector} in {wildcards.runid}"
    input:
        ideal_psl=lambda wc: patterns.output_ideal_psl_filename(
            config,
            hpge_detector=wc.hpge_detector,
            hpge_voltage=smk_load_hpge_cache()[wc.runid][wc.hpge_detector][
                "operational_voltage_in_V"
            ],
        ),
        electronics_model=lambda wc: patterns.output_elecmod_merged_filename(
            config, runid=wc.runid
        ),
    output:
        psl_file=patterns.output_realistic_psl_filename(config),
        plot_file=patterns.plot_realistic_psl_filename(config),
    log:
        patterns.log_realistic_psl_filename(config),
    benchmark:
        patterns.benchmark_realistic_psl_filename(config)
    script:
        "../src/legendsimflow/scripts/make_hpge_realistic_pulse_shape_lib.py"


rule merge_hpge_realistic_psls:
    """Merge HPGe realistic PSL in a single file.

    Copy the top-level LH5 objects from each individual detector realistic-PSL
    file into a single merged file using `h5copy`.

    Uses wildcard `runid`.
    """
    message:
        "Merging HPGe realistic psl files for {wildcards.runid}"
    input:
        lambda wc: aggregate.gen_list_of_realistic_psls(
            config, wc.runid, cache=smk_load_hpge_cache()
        ),
    output:
        patterns.output_realistic_psl_merged_filename(config),
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
    """Extract the HPGe current-pulse model.

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
    params:
        # track l200data so the rule reruns when it changes: inputs are
        # discovered dynamically and not listed above
        _l200data=config.paths.get("l200data", None),
    output:
        pars_file=temp(patterns.output_currmod_filename(config)),
        plot_file=patterns.plot_currmod_filename(config),
    log:
        patterns.log_currmod_filename(config),
    script:
        "../src/legendsimflow/scripts/extract_hpge_current_pulse_model.py"


rule merge_current_pulse_model_pars:
    """Merge the HPGe current-pulse model parameters in a single file per `runid`.

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


rule simulate_electron_gun:
    """Simulate mono-energetic electrons in the bulk of all HPGe detectors.

    The electron-gun simulations feed the A/E mean energy-dependence extraction
    (`extract_hpge_aoemean_energy_dependence`). This rule builds the production
    geometry from the experiment's template geometry configuration, renders the
    remage macro from the `aoemeanmod` macro template (GPS mono-energetic
    electrons confined in the bulk of all germanium sensitive volumes, see
    {func}`legendsimflow.commands.make_electron_gun_macro`), runs remage once
    per electron energy (900 to 2350 keV, the range corrected in data) through
    its Python
    API and plots the primary vertices as a check of the confinement.

    The `stp` output files (one per energy) and the geometry are temporary:
    they are deleted once the extraction has run for all runs.

    No wildcards are used.
    """
    message:
        "Simulating the electron gun"
    input:
        geom_config=patterns.geom_template_config_filename(config),
        template=electron_gun_macro_template(config),
    params:
        # rerun when the requested statistics change
        _primaries=electron_gun_primaries(config),
    output:
        geom=temp(patterns.output_electron_gun_geom_filename(config)),
        macro=patterns.output_electron_gun_macro_filename(config),
        stp_files=temp(aggregate.gen_list_of_electron_gun_stp_outputs(config)),
        plot_file=patterns.plot_electron_gun_vertices_filename(config),
    log:
        patterns.log_electron_gun_stp_filename(config),
    benchmark:
        patterns.benchmark_electron_gun_stp_filename(config)
    threads: 1
    script:
        "../src/legendsimflow/scripts/simulate_electron_gun.py"


def smk_extract_hpge_aoemean_energy_dependence_inputs(wildcards):
    """Inputs of the per-run A/E mean energy-dependence extraction.

    The electron-gun ``stp`` files (all energies and jobs) and the geometry they
    were simulated with, the HPGe modeling status cache and the PSD inputs of
    the run (drift-time maps, current-pulse model, realistic pulse-shape
    library), gated by the same hit-tier settings as in ``smk_hit_inputs``.
    """
    hit_settings = get_tier_settings(config, "hit")

    inputs = {
        "electron_stp_files": aggregate.gen_list_of_electron_gun_stp_outputs(config),
        "geom": patterns.output_electron_gun_geom_filename(config),
        "is_modelable": rules.cache_modelable_hpges.output.is_modelable,
    }
    if hit_settings.get("simulate_psd", True):
        inputs["hpge_dtmap"] = patterns.output_dtmap_merged_filename(
            config, runid=wildcards.runid
        )
        inputs["hpge_currmod"] = patterns.output_currmod_merged_filename(
            config, runid=wildcards.runid
        )
    if hit_settings.get("simulate_psd_with_psl", False):
        inputs["hpge_realistic_psl"] = patterns.output_realistic_psl_merged_filename(
            config, runid=wildcards.runid
        )
    return inputs


# realistic-PSL-based PSD simulation keeps a detector's PSL waveforms in
# memory, as in the hit tier (see hit.smk)
_aoemeanmod_resources = (
    {"mem_mb": 2300}
    if get_tier_settings(config, "hit").get("simulate_psd_with_psl", False)
    else {}
)


rule extract_hpge_aoemean_energy_dependence:
    """Extract the A/E mean energy dependence of the HPGe detectors of a run.

    Post-process the electron-gun simulations (`simulate_electron_gun`) for
    every modelable HPGe detector of the run with the same active-volume and
    PSD routines as `build_tier_hit` (drift-time maps, current-pulse model and
    realistic pulse-shape library of the run, no correction applied) and
    summarise the raw A/E distribution at each electron energy with its
    unbinned mode (the half-sample mode, insensitive to the low-side tail left
    by bremsstrahlung losses; its uncertainty comes from a bootstrap).
    The modes are fitted with a linear function of the energy, with the same
    model and procedure used on LEGEND-200 data, which is the A/E mean model
    applied by the `hit` tier.

    Two YAML files are produced. The model file, consumed by `build_tier_hit`,
    maps each detector to its `single_template` and `psl` models
    (`expression`, `pars`, `errs`), in the same format as the
    `aoemeanmod_default` hit-tier setting. The statistics file maps each
    detector to:

    | Key             | Type    | Description                                                                                                                 |
    | --------------- | ------- | --------------------------------------------------------------------------------------------------------------------------- |
    | `electron_gun`  | mapping | Per PSD method, arrays over the electron energies: `energy_in_keV`, `n_events`, `mode`, `mode_err`, `median`, `std`, `mad`. |
    | `mc_resolution` | mapping | Per PSD method, the standard deviation (in percent) of the simulated raw A/E at the electron energy closest to 2 MeV.      |

    A validation PDF shows, for each detector, the raw A/E distributions at
    each energy and the modes with the fitted model.

    Uses wildcard `runid`.
    """
    message:
        "Extracting the A/E mean energy dependence in {wildcards.runid}"
    input:
        unpack(smk_extract_hpge_aoemean_energy_dependence_inputs),
    output:
        pars_file=patterns.output_aoemeanmod_filename(config),
        stats_file=patterns.output_aoemeanmod_stats_filename(config),
        plot_file=patterns.plot_aoemeanmod_filename(config),
    log:
        patterns.log_aoemeanmod_filename(config),
    benchmark:
        patterns.benchmark_aoemeanmod_filename(config)
    resources:
        **_aoemeanmod_resources,
    script:
        "../src/legendsimflow/scripts/extract_hpge_aoemean_energy_dependence.py"


# whether superpulses are built one file per (run, detector) instead of a single
# file per detector accumulating all runs; drives both the producer rule below
# and the consumer `extract_electronics_model_pars`
_build_per_runid = get_par_settings(config, "superpulses").get("build_per_runid", False)


rule build_superpulses_from_data:
    """Build HPGe data superpulses (average waveforms per drift-time slice).

    Accumulate waveforms from LEGEND-200 data per detector (across all configured
    runs, or per run when `build_per_runid` is set) and average them within
    drift-time slices to form "super" pulses, after a self-similarity chi2 cut.

    The output LH5 file contains one group per drift-time slice, named
    `{detector}/dt_{lo}_{hi}_ns/`, each holding:

    | Field                  | Type     | Units             | Description                                                     |
    | ---------------------- | -------- | ----------------- | --------------------------------------------------------------- |
    | `charge_wf`            | `Array`  | ADC/cuspEmax      | Average normalised charge waveform.                             |
    | `current_wf`           | `Array`  | (ADC/cuspEmax)/ns | Average current waveform.                                       |
    | `charge_time_axis`     | `Array`  | ns                | Time axis for the charge waveform, aligned to `tp_aoe_max = 0`. |
    | `current_time_axis`    | `Array`  | ns                | Time axis for the current waveform.                             |
    | `drift_time_center`    | `Scalar` | ns                | Centre of the drift-time bin.                                   |
    | `drift_time_lo`        | `Scalar` | ns                | Lower bound of the drift-time bin.                              |
    | `drift_time_hi`        | `Scalar` | ns                | Upper bound of the drift-time bin.                              |
    | `energy_lo`            | `Scalar` | keV               | Lower bound of the energy selection.                            |
    | `energy_hi`            | `Scalar` | keV               | Upper bound of the energy selection.                            |
    | `detector`             | `Scalar` |                   | Detector name string.                                           |
    | `n_events_preliminary` | `Scalar` |                   | Waveforms before the chi2 cut.                                  |
    | `n_events_final`       | `Scalar` |                   | Waveforms after the chi2 cut.                                   |

    Uses wildcard `hpge_detector`, and `runid` if per-run superpulses are
    requested (`build_per_runid`).
    """
    message:
        "Building data superpulses for detector {wildcards.hpge_detector}"
    params:
        runids=lambda wc: (
            wc.runid
            if _build_per_runid
            else sorted(aggregate.gen_list_of_all_runids(config))
        ),
        # track l200data so the rule reruns when it changes: raw data files are
        # discovered dynamically and not listed as inputs
        _l200data=config.paths.get("l200data", None),
    output:
        superpulses=patterns.output_superpulses_filename(
            config, build_per_runid=_build_per_runid
        ),
        plots=patterns.plot_superpulses_filename(
            config, build_per_runid=_build_per_runid
        ),
    log:
        patterns.log_superpulses_filename(config, build_per_runid=_build_per_runid),
    script:
        "../src/legendsimflow/scripts/build_superpulses_from_data.py"


rule extract_electronics_model_pars:
    """Extract the HPGe electronics-response model.

    Fit the electronics transfer function (a Gaussian digitizer bandwidth `sigma`
    convolved with a causal exponential preamplifier decay `tau`) by comparing
    processed ideal PSL waveforms to the measured data superpulses. When a
    `default` key is present in the `elecmod` validity metadata the fit is
    bypassed and the metadata values are written directly.

    The (temporary, per-detector) output YAML, consumed by
    `merge_electronics_model_pars`, contains:

    | Key        | Type             | Description                                                                                                  |
    | ---------- | ---------------- | ----------------------------------------------------------------------------------------------------------- |
    | `detector` | str              | Detector name.                                                                                              |
    | `sigma`    | float            | Fitted Gaussian sigma of the digitizer bandwidth, in ns.                                                    |
    | `tau`      | float            | Fitted preamplifier decay constant, in ns.                                                                  |
    | `rms`      | float            | Best-fit mean RMS across slices.                                                                            |
    | `angle`    | str              | Crystal-axis angle tag used for PSL waveform selection (e.g. `"000"`).                                       |
    | `aoe_data` | float (optional) | Max current amplitude of the data superpulse in the highest drift-time slice; present with diagnostic plots. |
    | `aoe_mc`   | float (optional) | Max current amplitude of the simulation superpulse in the highest drift-time slice; present with plots.      |

    Only `sigma` and `tau` are required by the consumer
    (`convolve_hpge_ideal_pulse_shape_lib`).

    Uses wildcards `runid` and `hpge_detector`.
    """
    message:
        "Extracting electronics model for detector {wildcards.hpge_detector} in {wildcards.runid}"
    input:
        superpulses=lambda wc: (
            (
                patterns.output_superpulses_filename(
                    config, build_per_runid=True, runid=wc.runid
                )
                if _build_per_runid
                else patterns.output_superpulses_filename(config)
            )
            if patterns.compute_superpulses(config, runid=wc.runid)
            else []
        ),
        ideal_psl=lambda wc: patterns.output_ideal_psl_filename(
            config,
            hpge_detector=wc.hpge_detector,
            hpge_voltage=smk_load_hpge_cache()[wc.runid][wc.hpge_detector][
                "operational_voltage_in_V"
            ],
        ),
    output:
        pars_file=temp(patterns.output_elecmod_filename(config)),
        plot_file=patterns.plot_elecmod_filename(config),
        uniformity_plot_file=patterns.plot_superpulses_uniformity_filename(
            config, build_per_runid=True
        ),
    log:
        patterns.log_elecmod_filename(config),
    script:
        "../src/legendsimflow/scripts/extract_hpge_elec_response_model.py"


rule merge_electronics_model_pars:
    """Merge the HPGe electronics-response model parameters in a single file per `runid`.

    Collect the individual best-fit parameter files (one per detector) and
    write them into a single YAML file keyed by detector name.

    Uses wildcard `runid`.
    """
    message:
        "Merging electronics model parameters in {wildcards.runid}"
    input:
        lambda wc: aggregate.gen_list_of_elecmods(
            config, wc.runid, cache=smk_load_hpge_cache()
        ),
    params:
        # materialize the HPGe list here so the run: block can map file
        # indices back to detector names without calling back into the cache
        hpges=lambda wc: list(smk_load_hpge_cache()[wc.runid]),
    output:
        patterns.output_elecmod_merged_filename(config),
    run:
        import dbetto

        # NOTE: this is guaranteed to be sorted as in the input file list
        hpges = params.hpges

        out_dict = {}
        for i, f in enumerate(input):
            out_dict[hpges[i]] = dbetto.utils.load_dict(f)

        dbetto.utils.write_dict(out_dict, output[0])


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
    params:
        # track l200data so the rule reruns when it changes: the parameter
        # database is discovered dynamically and not listed above
        _l200data=config.paths.get("l200data", None),
    output:
        eresmod_file=patterns.output_eresmod_filename(config),
        aoeresmod_file=patterns.output_aoeresmod_filename(config),
        psdcuts_file=patterns.output_psdcuts_filename(config),
    log:
        patterns.log_eresmod_filename(config),
    script:
        "../src/legendsimflow/scripts/pars/extract_hpge_observables_models.py"
