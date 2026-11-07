# Output files

This page describes the files a production writes and what they contain. The
steps that write them are described in [](pipeline.md), and the output folders
(`paths.tier`, `paths.pars`, ...) are set in the configuration file (see
[](setup.md)). The simulation and event files use the
[LEGEND HDF5 (LH5) format](https://legend-exp.github.io/legend-data-format-specs/dev/hdf5/);
the parameters of the `par` step are YAML or LH5 files.

(file-naming)=

## File names

Each step writes to its own folder, `paths.tier.<tier>`:

| Step                              | File                                                                    |
| --------------------------------- | ----------------------------------------------------------------------- |
| `vtx`, `stp`, `opt`, `hit`, `evt` | `{simid}/{experiment}-{simid}-job_{jobid}-tier_{tier}.lh5`, one per job |
| `cvt`                             | `{experiment}-{simid}-tier_cvt.lh5`, one per simulation                 |
| `pdf`                             | `{experiment}-{simid}-tier_pdf.lh5`, one per simulation                 |

`experiment` is the experimental configuration, `simid` the simulation (see
{ref}`simconfig.yaml`) and `jobid` the zero-padded job number. Hyphens separate
these fields, which is why simids should use underscores. The GDML geometry of a
simulation is `paths.geom/{experiment}-{simid}-tier_stp-geom.gdml` (see
{ref}`geom-config.yaml`). The `par` step writes its files below `paths.pars`
(see {ref}`par-output-files`). Control plots go to `paths.tier.<tier>/plots/`,
and logs and benchmarks to `paths.log` and `paths.benchmarks`.

(par-output-files)=

## `par` step

The parameters of the `par` step (see [](par.md)), below `paths.pars` unless
noted, per run (`{runid}`) and detector (`{detector}`):

| Parameter                  | File                                                                                                          |
| -------------------------- | ------------------------------------------------------------------------------------------------------------- |
| run partitions             | `simstat/partitions_{simid}.yaml`, see {ref}`run-partitioning`                                                |
| energy resolution          | `hpge/eresmod/{runid}-model.yaml`, see {ref}`hpge-eresmod-extraction`                                         |
| A/E resolution             | `hpge/aoeresmod/{runid}-model.yaml`, see {ref}`hpge-aoeresmod-extraction`                                     |
| PSD cuts                   | `hpge/psdcuts/{runid}-psd-cuts.yaml`, see {ref}`hpge-psdcuts-extraction`                                      |
| current-pulse model        | `hpge/currmod/{runid}-model.yaml`, see {ref}`hpge-currmod-extraction`                                         |
| A/E mean energy dependence | `hpge/aoemeanmod/model.yaml`, for all runs, see {ref}`hpge-aoemeanmod-extraction`                             |
| drift-time maps            | `paths.dtmaps/{runid}-hpge-drift-time-maps.lh5`, see {ref}`hpge-dtmap-extraction`                             |
| ideal pulse shapes         | `hpge/psl/ideal/singles/{detector}-{voltage}V-hpge-pulse-shape-lib.lh5`, see {ref}`hpge-ideal-psl-extraction` |
| data superpulses           | `hpge/superpulses/{detector}-superpulses.lh5`, see {ref}`hpge-superpulses-extraction`                         |
| electronics response       | `hpge/elecmod/{runid}-model.yaml`, see {ref}`hpge-elecmod-extraction`                                         |
| realistic pulse shapes     | `hpge/psl/realistic/{runid}-hpge-pulse-shape-lib.lh5`, see {ref}`hpge-realistic-psl-extraction`               |
| tuned impurity profiles    | `hpge/impurity.yaml`, see [](impurity-tuning.md)                                                              |

Each file holds all the detectors of the run; the files of single detectors that
they combine sit in `singles/` subfolders next to them. With `build_per_runid`
(see {ref}`superpulses-settings-meta`), the superpulses are written per run, as
`{runid}-{detector}-superpulses.lh5`. Which files a production builds depends on
its settings: the pulse-shape files only with PSD simulation, the A/E mean only
with the two-pass correction, the tuned profiles only in an impurity tuning (see
{ref}`hit-tier-settings`).

(par-detinfo)=

## Detector information

The detector information that the steps need per run is read from the LEGEND
metadata once and stored in `paths.pars/detinfo/`, one YAML file per quantity,
each a mapping `runid -> detector -> value`:

```yaml
# usability.yaml
l200-p03-r000-phy:
  V00048A: on
  B00035B: ac
  S002: off
```

| File                              | Detectors | Value                         | Meaning                                                                           |
| --------------------------------- | --------- | ----------------------------- | --------------------------------------------------------------------------------- |
| `usability.yaml`                  | all       | `on`, `ac`, `off`             | the status of the detector in the run                                             |
| `rawid.yaml`                      | all       | integer                       | the channel number in the data; also lists the runs random coincidences come from |
| `psd_usability.yaml`              | HPGe      | `valid`, `present`, `missing` | the status of the PSD of the detector in the data; `valid` when not set           |
| `crystal_metadata_usability.yaml` | HPGe      | status string or `null`       | whether the crystal metadata of the detector is good enough for the simulation    |
| `is_modelable.yaml`               | HPGe      | `true`, `false`               | whether the detector gets a pulse-shape model (see {ref}`hpge-modeling-criteria`) |
| `operational_voltage_in_V.yaml`   | HPGe      | integer or `null`             | the operating voltage of the detector in the run; `null` without one              |
| `hpge_ssd_modeling.yaml`          | HPGe      | mapping                       | the field simulation of each modelled detector (see below)                        |

`hpge_ssd_modeling.yaml` records, for each modelled detector, the factor that
scaled its impurity density (`impurity_scaling_factor`), the target depletion
voltage (`measured_depletion_voltage_in_V`) and the simulated one before and
after the scaling (`simulated_depletion_voltage_raw_in_V`,
`simulated_depletion_voltage_in_V`), all in V (see {ref}`hpge-psl-overview`).

The files are not rebuilt when the metadata changes (see [](meta.md)).

(hit-tier-output)=

## `hit` step

The `hit` file of a job (see [](hit.md)) holds one table per HPGe detector,
`/hit/{detector}`, with one row per hit, and a TCM (see {ref}`tcm-output`). The
events of each run partition (see {ref}`run-partitioning`) follow each other in
the table.

| Field           | Type    | Units | Description                                                                                             |
| --------------- | ------- | ----- | ------------------------------------------------------------------------------------------------------- |
| `evtid`         | integer |       | the identifier of the simulated primary                                                                 |
| `t0`            | float   | ns    | the time of the hit, from the start of the primary                                                      |
| `energy`        | float   | keV   | the energy, after the active volume model and the energy resolution                                     |
| `period`, `run` | integer |       | the period and run of the partition, e.g. 3 and 1 for `l200-p03-r001-phy`                               |
| `usability`     | integer |       | the status of the detector in the run: 0 on, 1 ac, 2 off                                                |
| `psd_usability` | integer |       | the status of its PSD in the data: 0 valid, 1 present, 2 missing                                        |
| `is_valid_sim`  | bool    |       | whether the crystal metadata of the detector is good enough for the simulation (see {ref}`par-detinfo`) |

A detector missing from the detector information of the run counts as on, with
valid PSD, and `is_valid_sim` false. With PSD simulation (see
{ref}`hit-tier-settings`), the `psd` subtable has one subtable per method,
`single_temp` and `pulse_lib` (see [](hit.md)):

| Field             | Type  | Units | Method      | Description                                                                         |
| ----------------- | ----- | ----- | ----------- | ----------------------------------------------------------------------------------- |
| `aoe_raw`         | float |       | both        | the maximum of the current divided by the energy                                    |
| `aoe_corr`        | float |       | both        | the A/E corrected for its energy dependence (see {ref}`hpge-aoemeanmod-extraction`) |
| `aoe`             | float |       | both        | the A/E classifier, `(aoe_corr - 1) / A/E resolution`                               |
| `drift_time_amax` | float | ns    | both        | the drift time at the maximum of the current                                        |
| `is_single_site`  | bool  |       | both        | the classifier is above the low PSD cut                                             |
| `is_high_aoe`     | bool  |       | `pulse_lib` | the classifier is above the high PSD cut                                            |
| `is_bb_like`      | bool  |       | `pulse_lib` | single-site and not high A/E                                                        |

A detector without a PSD model (no drift-time map or current-pulse template for
the single template, no pulse-shape library for the other) has `NaN` A/E fields
and false flags (see {ref}`build-tier-hit-hpge`).

With `store_surface_steps` (see {ref}`hit-tier-settings`), the `surface`
subtable holds the energy depositions near the n+ surface, as deposited, before
the active volume model:

| Field               | Type           | Units | Description                                                                 |
| ------------------- | -------------- | ----- | --------------------------------------------------------------------------- |
| `edep`              | list of floats | keV   | the energy of each step closer than `surface_depth_in_mm` to the n+ surface |
| `distance_to_nplus` | list of floats | mm    | the distance of each of these steps to the n+ surface                       |
| `energy_bulk`       | float          | keV   | the summed energy of all the other steps                                    |

(opt-tier)=

## `opt` step

The `opt` file of a job (see [](opt.md)) holds the photoelectrons of the SiPMs
below `/hit/`: one table per SiPM with per-SiPM optical maps or tracked photons,
or one table `/hit/spms` for all SiPMs together with a single map, and a TCM
(see {ref}`tcm-output`). With optical maps, a row is an event window with light
from a light source; with tracked photons, a hit of the SiPM.

| Field           | Type           | Units          | Description                                                                                                    |
| --------------- | -------------- | -------------- | -------------------------------------------------------------------------------------------------------------- |
| `evtid`         | integer        |                | the identifier of the simulated primary                                                                        |
| `t0`            | float          | ns             | the time of the row: the earliest hit of the window, or the SiPM hit                                           |
| `dt`            | list of floats | ns             | the times of the photoelectrons, from `t0`, after merging the ones closer than the time resolution             |
| `energy`        | list of floats | photoelectrons | their amplitudes, smeared with the single-photoelectron resolution                                             |
| `is_saturated`  | bool           |                | the row reached the cap on photoelectrons (see {ref}`opt-tier-settings`)                                       |
| `expected_pes`  | float          |                | the expected number of photoelectrons at unit SiPM efficiency, with optical maps and `store_expected_pes` only |
| `period`, `run` | integer        |                | the period and run of the partition                                                                            |
| `usability`     | integer        |                | the status of the SiPM in the run (0 on, 1 ac, 2 off); always on in `/hit/spms`                                |

(evt-tier)=

## `evt` step

The `evt` file of a job (see [](evt.md)) holds one row per event in the table
`/evt`, organized like the `evt` tier of the LEGEND-200 data, and at its root:

- `tcm`, the TCM of the events (see {ref}`tcm-output`);
- `detector_uids`, the identifier of each detector by name, to read the `rawid`
  fields below;
- `number_of_simulated_events`, the number of simulated primaries of the job.

Without HPGe or SiPM data (`skip_hit`, `skip_opt`, see
{ref}`evt-tier-settings-meta`), the `geds` or `spms` subtable and the
corresponding coincidence flag are absent.

### `trigger/`

| Field                    | Type  | Units | Description                                                                              |
| ------------------------ | ----- | ----- | ---------------------------------------------------------------------------------------- |
| `evtid`, `period`, `run` | float |       | the primary, period and run of the first hit of the event                                |
| `timestamp`              | float | ns    | the time of the first hit of the event, from the start of the primary (double precision) |

The first hit is taken from the `hit` file, or from the `opt` file without HPGe
data. An event without any `hit`-file hit has `NaN` in these fields.

### `geds/`

The HPGe hits above `geds_energy_thr_kev` in detectors that are not off:

| Field                             | Type             | Units | Description                                                                                                   |
| --------------------------------- | ---------------- | ----- | ------------------------------------------------------------------------------------------------------------- |
| `energy`                          | list of floats   | keV   | the energy of each hit                                                                                        |
| `rawid`                           | list of integers |       | the detector of each hit, by its identifier in the simulated geometry (see `detector_uids`)                   |
| `hit_idx`                         | list of integers |       | the row of each hit in the `hit` file of the job                                                              |
| `multiplicity`                    | integer          |       | the number of hits, AC detectors included                                                                     |
| `energy_sum`                      | float            | keV   | the summed energy of the hits in detectors that are on                                                        |
| `quality/is_good_channel`         | list of bools    |       | the detector of each hit is on (not AC)                                                                       |
| `psd/is_good`                     | list of bools    |       | the PSD of the detector is valid in the data                                                                  |
| `psd/is_valid_sim`                | list of bools    |       | the crystal metadata is good enough for the simulation                                                        |
| `psd/{single_temp,pulse_lib}/...` | lists            |       | the PSD fields of the `hit` step (see {ref}`hit-tier-output`), plus `has_aoe`, true when the A/E is not `NaN` |
| `surface/...`                     | lists            |       | the surface depositions of the `hit` step, one list per hit, with `store_surface_steps` only (see [](evt.md)) |

### `spms/`

Every SiPM channel that is not off in the run, in ascending identifier order,
whether it has light or not:

| Field                  | Type                    | Units              | Description                                                                                  |
| ---------------------- | ----------------------- | ------------------ | -------------------------------------------------------------------------------------------- |
| `rawid`                | list of integers        |                    | the channels, by their identifier in the simulated geometry                                  |
| `energy`               | list of lists of floats | photoelectrons     | the amplitudes of the pulses of each channel above `spms_energy_thr_pe`                      |
| `time`                 | list of lists of floats | ns                 | their times, from `trigger/timestamp`                                                        |
| `is_saturated`         | list of bools           |                    | a row of the channel reached the cap on photoelectrons of the `opt` step                     |
| `hit_idx`              | list of integers        |                    | the first row of the channel in the `opt` file of the job, or -1                             |
| `expected_pes`         | list of floats          |                    | the expected photoelectrons of each channel, with optical maps and `store_expected_pes` only |
| `multiplicity`         | integer                 |                    | the number of channels with at least one pulse, random coincidences included                 |
| `energy_sum`           | float                   | photoelectrons     | the summed amplitude of all pulses, random coincidences included                             |
| `rc_energy`, `rc_time` | list of lists of floats | photoelectrons, ns | the random-coincidence pulses, with `add_random_coincidences` only; no threshold is applied  |

### `coincident/`

| Field  | Type | Description                                                                                                    |
| ------ | ---- | -------------------------------------------------------------------------------------------------------------- |
| `geds` | bool | the event has at least one HPGe hit                                                                            |
| `spms` | bool | the LAr veto: the SiPM multiplicity or amplitude sum reaches its threshold (see {ref}`evt-tier-settings-meta`) |

(cvt-tier)=

## `cvt` step

The `cvt` file of a simulation (see [](evt.md)) concatenates the `evt` and `tcm`
tables of all the jobs, in job order. `number_of_simulated_events` is the sum
over the jobs, and `detector_uids` the union. The `hit_idx` fields and the
`row_in_table` of the TCM keep the values of their job, so they do not point to
valid rows of the concatenated file.

(tcm-output)=

## Time-coincidence map (TCM)

The TCM groups the hits of the detector tables into events (see
{ref}`detector-tables-and-events`). The `stp` file written by _remage_ has one,
and the `opt` and `hit` files rebuild it for their own tables: hits of the same
primary are in the same event when each is less than 10 µs after the previous
one. The `evt` file holds the combined TCM of the `stp` and `opt` files, whose
entries point to the tables of the `hit` and `opt` files of the job. It has one
row per row of the `stp` TCM.

| Field          | Type             | Description                                                    |
| -------------- | ---------------- | -------------------------------------------------------------- |
| `table_key`    | list of integers | the detector table of each hit of the event, by its identifier |
| `row_in_table` | list of integers | the row of each hit in that table                              |

(pdf-tier)=

## `pdf` step

The `pdf` file of a simulation (see [](pdf.md)) holds the histograms of the HPGe
energy below `pdf/`, and `nr_sim_events`, the number of simulated primaries of
the simulation, used to normalise them. The energy histograms of the
multiplicity-1 events are at

```text
pdf/mul1/<lar_category>/<aoe_category>[/<fccd_category>]/<group>
```

with 1 keV bins from 0 to 6000 keV, one per detector group (`all` and the groups
of {ref}`pdf-tier-settings`). The two-dimensional histogram of the
multiplicity-2 events is `pdf/mul2/all`, with the lower energy on the first
axis. The categories are:

| Level | Categories                                                                            |
| ----- | ------------------------------------------------------------------------------------- |
| LAr   | `is_spms_coincident`, `not_spms_coincident`                                           |
| A/E   | `is_high_aoe` (pulse-shape library only), `is_ss_aoe`, `is_low_aoe`, `is_invalid_aoe` |
| FCCD  | `is_fccd_tagged`, `not_fccd_tagged` (optional)                                        |

Their definitions are on [](pdf.md) and in {mod}`legendsimflow.selections` (see
[](analysis.md)). The categories of a level do not overlap, so any other
selection is a sum of them.

Without HPGe data, the file holds only `nr_sim_events`.

(pdf-tier-data)=

### LEGEND-200 data

The histograms of the LEGEND-200 data (see [](pdf.md)) are
`paths.tier.pdf/l200data/{runid}-tier_pdf.lh5`, one per run, and their sum
`paths.tier.pdf/l200data/{experiment}-l200data-tier_pdf.lh5`. They have the
layout of the simulated ones, without the FCCD level, and store `livetime_in_s`,
the livetime of the runs in s, instead of `nr_sim_events`.
