# Output tier fields

This page documents the fields (columns) produced by each post-processing tier
of the Simflow. All output files use the
[LEGEND HDF5 (LH5) format](https://legend-exp.github.io/legend-data-format-specs/dev/hdf5/).

(par-output-files)=

## `par` tier — HPGe PSD parameter files

The drift-time-map and pulse-shape-library rules of the `par` tier write the
per-detector / per-run parameter files consumed by `build_tier_hit` for PSD
simulation. See [](pipeline.md) for how each is produced, and the rule reference
linked there for the internal field structure. Their on-disk locations are:

| Output                              | Location                                                                                    |
| ----------------------------------- | ------------------------------------------------------------------------------------------- |
| Drift-time map (merged per run)     | `{config.paths.dtmaps}/{runid}-hpge-drift-time-maps.lh5`                                    |
| Ideal PSL (per detector/voltage)    | `{config.paths.pars}/hpge/psl/ideal/singles/{detector}-{voltage}V-hpge-pulse-shape-lib.lh5` |
| Data superpulses                    | `{config.paths.pars}/hpge/superpulses/{detector}-superpulses.lh5`                           |
| Electronics-response model (merged) | `{config.paths.pars}/hpge/elecmod/{runid}-model.yaml`                                       |
| Realistic PSL (merged per run)      | `{config.paths.pars}/hpge/psl/realistic/{runid}-hpge-pulse-shape-lib.lh5`                   |

The data-superpulse layout switches to one file per `(runid, detector)` pair
(`{runid}-{detector}-superpulses.lh5`) when `build_per_runid` is set (see
{ref}`superpulses-settings-meta`). The drift-time map and realistic PSL also
have per-detector `singles/` files that the merge step combines.

(par-detinfo)=

## `par` tier — detector-info cache (`pars/detinfo/`)

Post-processing needs per-detector, per-run information queried from
`legend-metadata` (channel-map status, diode and crystal records). Those lookups
are slow, so the `par` tier caches the results once as a set of small YAML files
under `{config.paths.pars}/detinfo/`, one file per _flag_, and the `opt`, `hit`,
and modeling rules read them instead of re-querying the metadata.

Every file is named `{flag}.yaml` and holds a two-level mapping
`runid -> detector -> value`:

```yaml
# usability.yaml
l200-p03-r000-phy:
  V00048A: on
  B00035B: ac
  S002: off
l200-p03-r001-phy:
  V00048A: on
  ...
```

A detector appears under a flag only when that flag applies to it: SiPM channels
carry `usability` and `rawid` alone, while `psd_usability`,
`crystal_metadata_usability`, `is_modelable`, and `operational_voltage_in_V` are
germanium-only.

The [`cache_detector_usabilities`](../api/snakemake_rules.md) rule writes the
data-quality flags for **all** deployed detectors:

| File                              | Applies to  | Value                | Description                                                                                                                                                             |
| --------------------------------- | ----------- | -------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `usability.yaml`                  | geds + SiPM | `on`, `off`, `ac`, … | Analysis usability from the channel-map status (`analysis.usability`).                                                                                                  |
| `rawid.yaml`                      | geds + SiPM | integer              | Rawid from the channel map (`daq.rawid`). Changes when channels are recabled. Also lists the runs that random coincidences are drawn from (`random_coincidence_runid`). |
| `psd_usability.yaml`              | geds        | `valid`, …           | PSD usability from `analysis.psd.status.low_aoe`; defaults to `valid` when the field is absent.                                                                         |
| `crystal_metadata_usability.yaml` | geds        | `valid`, …, `null`   | Usability of the crystal metadata required for modeling (the crystal-slice `status`). `null` when the information is unavailable in the metadata.                       |

The [`cache_modelable_hpges`](../api/snakemake_rules.md) checkpoint writes the
modeling flags for **every deployed germanium** detector. Every run of the
Simflow is listed, with no detector if no HPGe is deployed in it. The checkpoint
is not run when the `evt` tier skips the `hit` tier (`skip_hit`):

| File                            | Value            | Description                                                                                                                                     |
| ------------------------------- | ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `is_modelable.yaml`             | `true` / `false` | Whether the detector is suitable for drift-time-map and current-pulse modeling. See {ref}`hpge-modeling-criteria` for the eligibility criteria. |
| `operational_voltage_in_V.yaml` | integer / `null` | Operational bias voltage read from the parameters database. `null` for detectors with no bias (e.g. off detectors).                             |

Unlike the flags above (queried from `legend-metadata`), the
[`aggregate_hpge_ssd_modeling_info`](../api/snakemake_rules.md) rule writes a
single file whose values are **outputs of the SSD simulation** run for each
drift-time map. Its `runid -> detector` leaves are a nested mapping rather than
a scalar:

| File                     | Applies to | Description                                                                         |
| ------------------------ | ---------- | ----------------------------------------------------------------------------------- |
| `hpge_ssd_modeling.yaml` | geds       | SSD-modeling provenance for each modelable detector: the four scalars listed below. |

| Key                                    | Type           | Units | Description                                                                                                                             |
| -------------------------------------- | -------------- | ----- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `impurity_scaling_factor`              | float / `null` | —     | Dimensionless factor applied to the crystal impurities to match the measured depletion voltage. `null` when no rescaling was performed. |
| `measured_depletion_voltage_in_V`      | float / `null` | V     | Measured depletion voltage from `legend-metadata` (`characterization.l200_site.depletion_voltage_in_V`). `null` when absent.            |
| `simulated_depletion_voltage_raw_in_V` | float          | V     | Depletion voltage estimated from the simulation before the impurity rescaling.                                                          |
| `simulated_depletion_voltage_in_V`     | float          | V     | Depletion voltage estimated from the simulation after the impurity rescaling and bias adjustment.                                       |

This file is produced together with the drift-time maps, so it is written only
when the `hit` tier is built with PSD enabled (`simulate_psd`, see
{ref}`hit-tier-settings`).

(hit-tier-output)=

## `hit` tier — HPGe detector post-processing

The `hit` tier applies detector response models (energy resolution, pulse-shape
discrimination) to the raw
[`stp`-tier](https://remage.readthedocs.io/en/stable/manual/output.html)
simulation output. Each HPGe detector is processed independently; output tables
are stored under `/hit/{detector_name}/` in the LH5 file. Each row corresponds
to a single _remage_ hit in one detector — use the `evtid` column and the
[`evt` tier](evt-tier) to group hits into physics events.

### Inherited fields

These fields are carried over from the
[`stp` tier](https://remage.readthedocs.io/en/stable/manual/output.html):

| Field   | Type    | Units | Description                                                                                    |
| ------- | ------- | ----- | ---------------------------------------------------------------------------------------------- |
| `evtid` | `Array` | —     | Event identifier, shared across all detectors hit in the same event.                           |
| `t0`    | `Array` | ns    | Time of the first energy deposition in the detector relative to the start of the Geant4 event. |

### Added fields

| Field           | Type    | Units | Description                                                                                                                                                                                                                                                                                     |
| --------------- | ------- | ----- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `energy`        | `Array` | keV   | Reconstructed energy after smearing with the detector energy resolution. Computed from the sum of active energy depositions (weighted by the dead-layer activeness model).                                                                                                                      |
| `period`        | `Array` | —     | Data-taking period number extracted from the run identifier (numeric encoding).                                                                                                                                                                                                                 |
| `run`           | `Array` | —     | Data-taking run number extracted from the run identifier (numeric encoding).                                                                                                                                                                                                                    |
| `usability`     | `Array` | —     | Encoded detector usability status for this run (e.g. `on`, `off`, `ac`). Decode with {func}`legendsimflow.metadata.decode_usability`. See the detector status flags in `legend-metadata/datasets/statuses`.                                                                                     |
| `psd_usability` | `Array` | —     | Encoded PSD usability flag (e.g. `valid`). Indicates whether PSD parameters are valid in LEGEND-200 data for this detector and run. Decode with {func}`legendsimflow.metadata.decode_psd_usability`.                                                                                            |
| `is_valid_sim`  | `Array` | —     | Boolean. `True` when the crystal metadata needed to model this detector is usable (`crystal_metadata_usability` is `valid`, see {ref}`par-detinfo`), so the simulated detector response can be trusted. `False` otherwise. Reshaped into `geds/psd/is_valid_sim` by the [`evt` tier](evt-tier). |

The `hit` tier can add HPGe pulse-shape-discrimination (PSD) fields in the `psd`
subtable, which groups one method-specific subtable per simulation method,
selected by the metadata settings in {ref}`hit-tier-settings`:

- `psd/single_temp` — single-template A/E simulation, present only when
  `simulate_psd: True` (the default).
- `psd/pulse_lib` — pulse-shape-library (PSL) based simulation, present only
  when `simulate_psd_with_psl: True` (see {ref}`hpge-psl-overview`).

:::{note}

The two subtables hold the same-named A/E fields, but their **meaning differs**:
the `single_temp` values come from a single per-detector current-pulse template,
while the `pulse_lib` values come from the per-pixel pulse-shape library (see
{ref}`hpge-psl-overview`).

:::

| Field             | Type    | Units | Subtable                   | Description                                                                                                                                                                                                                                                                                            |
| ----------------- | ------- | ----- | -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `drift_time_amax` | `Array` | ns    | `single_temp`, `pulse_lib` | Drift time at the maximum-current (A) position of the simulated current pulse. Set to `NaN` when no drift-time map or current-pulse model is available for the detector.                                                                                                                               |
| `aoe_raw`         | `Array` |       | `single_temp`, `pulse_lib` | Raw A/E value: maximum current amplitude divided by energy. The maximum current (A) is obtained from the simulated current pulse and includes electronic noise effects.                                                                                                                                |
| `aoe_corr`        | `Array` |       | `single_temp`, `pulse_lib` | Energy-corrected A/E value, obtained by correcting `aoe_raw` for the observed energy dependence.                                                                                                                                                                                                       |
| `aoe`             | `Array` |       | `single_temp`, `pulse_lib` | A/E classifier value: `(aoe_corr - 1) / aoe_resolution`. Used for pulse-shape discrimination (PSD). Set to `NaN` when PSD simulation is not available (i.e., when the drift-time map or current-pulse model is missing). This is distinct from the usability flags that track LEGEND-200 data quality. |
| `is_single_site`  | `Array` |       | `single_temp`, `pulse_lib` | Boolean PSD flag. `True` when the A/E classifier `aoe` exceeds the lower single-site cut `psdcuts.aoe.low_side` (extracted from LEGEND-200 data).                                                                                                                                                      |
| `is_bb_like`      | `Array` |       | `pulse_lib`                | Boolean PSD flag. `True` for $0\nu\beta\beta$-like single-site events: `aoe` above `psdcuts.aoe.low_side` and not above the upper cut `psdcuts.aoe.high_side`.                                                                                                                                         |
| `is_high_aoe`     | `Array` |       | `pulse_lib`                | Boolean PSD flag. `True` when `aoe` exceeds the upper cut `psdcuts.aoe.high_side` (high-A/E events, e.g. surface or $\alpha$).                                                                                                                                                                         |

With `store_surface_steps: true` (see {ref}`hit-tier-settings`), the `surface`
subtable holds the energy depositions (steps) closer than `surface_depth_in_mm`
to the n+ surface, and the summed energy of all other steps. The energies are as
deposited in the germanium, before any dead-layer or resolution model.

| Field               | Type              | Units | Description                                                                     |
| ------------------- | ----------------- | ----- | ------------------------------------------------------------------------------- |
| `edep`              | `VectorOfVectors` | keV   | Energy of each step closer than `surface_depth_in_mm` to the n+ surface.        |
| `distance_to_nplus` | `VectorOfVectors` | mm    | Distance of each of these steps to the n+ surface, in the same order as `edep`. |
| `energy_bulk`       | `Array`           | keV   | Summed energy of all other steps of the hit.                                    |

(opt-tier)=

## `opt` tier — optical (SiPM) post-processing

The `opt` tier is at the same conceptual level as the `hit` tier: it performs
detector-wise post-processing, but for SiPMs instead of HPGe detectors. It
applies the optical map convolution and photoelectron (PE) response models to
the scintillator output from the `stp` tier. When using per-SiPM optical maps, a
separate table is written for each SiPM channel under `/hit/{sipm_name}/`; when
using a single summed map (the default), all SiPM channels are aggregated into a
single `/hit/spms/` table. Each row corresponds to an `stp`-tier hit entry
(identified by `evtid`) — not to a single physics event: with the optical map,
an `stp` event window holding hits in one or more scintillator tables (liquid
argon, PEN, ...), whose light is pooled; with `light_source: tracked_photons`, a
hit of the SiPM itself (see {ref}`opt-tier-settings`). In the latter case every
SiPM of the geometry has a table, empty if no photon reached it.

### Inherited fields

| Field   | Type    | Units | Description                                                                                                                                      |
| ------- | ------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `evtid` | `Array` | —     | Event identifier, shared across all detectors hit in the same event.                                                                             |
| `t0`    | `Array` | ns    | Time of the first scintillator hit in the window (first photon in the SiPM hit with tracked photons), relative to the start of the Geant4 event. |

### Added fields

| Field          | Type              | Units | Description                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| -------------- | ----------------- | ----- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `dt`           | `VectorOfVectors` | ns    | Photoelectron arrival times relative to the `t0` of this row, after resolution smearing and photoelectron clustering to simulate the timing resolution of the SiPM. With the optical map or a light model, a photoelectron arrives at the time of its energy deposition plus the scintillation delay. With tracked photons, the photon transport and wavelength shifting add to this, and the first `dt` of each row is 0. Variable-length per row. |
| `energy`       | `VectorOfVectors` | —     | Photoelectron amplitudes (relative units), after PE resolution smearing. Variable-length array matching `dt`.                                                                                                                                                                                                                                                                                                                                       |
| `is_saturated` | `Array`           | —     | Boolean flag. `True` when the number of detected photoelectrons exceeds a maximum PE-per-hit cap, indicating SiPM saturation.                                                                                                                                                                                                                                                                                                                       |
| `expected_pes` | `Array`           | —     | _(optional)_ Expected number of photoelectrons per row at unit channel efficiency, before the PE-per-hit cap. Present only when `store_expected_pes` is enabled, with the optical map or a light model.                                                                                                                                                                                                                                             |
| `period`       | `Array`           | —     | Data-taking period number extracted from the run identifier (numeric encoding).                                                                                                                                                                                                                                                                                                                                                                     |
| `run`          | `Array`           | —     | Data-taking run number extracted from the run identifier (numeric encoding).                                                                                                                                                                                                                                                                                                                                                                        |
| `usability`    | `Array`           | —     | Encoded SiPM channel usability status for this run. Decode with {func}`legendsimflow.metadata.decode_usability`.                                                                                                                                                                                                                                                                                                                                    |

(evt-tier)=

## `evt` tier — event-level output

The `evt` tier merges HPGe (`hit`) and SiPM (`opt`) data into a unified
event-level structure. It uses the time-coincidence map (TCM) to associate hits
across detector subsystems into physics events. The output table is stored under
`/evt/` in the LH5 file and is organized into subtables. The structure is
designed to mirror the `evt` tier of the actual LEGEND-200 data (produced by
[pygama](https://legend-pydataobj.readthedocs.io)) as closely as possible.

:::{note}

The `geds/` and `coincident/geds` subtables are absent when `skip_hit: true` is
set in the evt tier settings (HPGe tier skipped). Similarly, the `spms/` and
`coincident/spms` subtables are absent when `skip_opt: true` is set (SiPM/LAr
tier skipped). See {ref}`evt-tier-settings-meta` for details.

:::

Each `evt` file also carries a root-level `number_of_simulated_events` scalar,
forwarded verbatim from the `number_of_simulated_events` field that _remage_
writes into the `stp` file. The `cvt` tier sums these per-job counts into a
single `number_of_simulated_events` scalar, which the `pdf` tier reads back as
`nr_sim_events`.

### `trigger/` — event metadata

Constant fields identifying each event.

| Field       | Type    | Units | Description                                                                                                                                                                                                                                                                                                                                                                                 |
| ----------- | ------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `evtid`     | `Array` | —     | Event identifier.                                                                                                                                                                                                                                                                                                                                                                           |
| `period`    | `Array` | —     | Data-taking period number.                                                                                                                                                                                                                                                                                                                                                                  |
| `run`       | `Array` | —     | Data-taking run number.                                                                                                                                                                                                                                                                                                                                                                     |
| `timestamp` | `Array` | ns    | `t0` of the first HPGe hit, relative to the start of the simulated Geant4 event. When the `hit` tier is skipped (`skip_hit: true`), `t0` of the first `opt`-tier hit: a scintillator hit window with the optical map, a SiPM hit with tracked photons. Used as the event's timestamp. Stored in double precision, because decays late in a chain happen hours after the start of the event. |

### `geds/` — HPGe detector array

Per-event arrays collecting HPGe hits that pass the energy threshold (25 keV)
and are from non-OFF detectors.

| Field          | Type              | Units | Description                                                                                                                                                                                               |
| -------------- | ----------------- | ----- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `energy`       | `VectorOfVectors` | keV   | Hit energies from ON and AC detectors above threshold. Variable-length per event.                                                                                                                         |
| `energy_sum`   | `Array`           | keV   | Summed energy from ON detectors only (excludes AC). Scalar per event.                                                                                                                                     |
| `rawid`        | `VectorOfVectors` | —     | Detector UID in the simulated geometry for each hit. Named `rawid` to match LEGEND-200 data; equal to the data rawid only if the geometry was built from the same channel map. Variable-length per event. |
| `hit_idx`      | `VectorOfVectors` | —     | Row index in the `hit`-tier table, for looking up additional hit-level fields. Variable-length per event.                                                                                                 |
| `multiplicity` | `Array`           | —     | Number of HPGe hits above threshold per event. Scalar per event.                                                                                                                                          |

#### `geds/quality/` — data-quality flags

| Field             | Type              | Description                                                                                 |
| ----------------- | ----------------- | ------------------------------------------------------------------------------------------- |
| `is_good_channel` | `VectorOfVectors` | Boolean. `True` if the detector usability is ON (not AC or OFF). Variable-length per event. |

#### `geds/psd/` — PSD fields

The optional `geds/psd` subtable holds HPGe pulse-shape-discrimination fields.
Flags common to every PSD method live directly under `geds/psd`, while
method-specific fields are grouped in dedicated subtables, one per simulation
method:

- `geds/psd/single_temp` — single-template A/E simulation, present when
  `simulate_psd: True` in {ref}`hit-tier-settings`.
- `geds/psd/pulse_lib` — pulse-shape-library (PSL) based simulation, present
  when `simulate_psd_with_psl: True` (see {ref}`hpge-psl-overview`).

The `geds/psd` subtable itself is present whenever at least one of the two
methods is enabled. The two method subtables carry the same-named A/E fields but
their **meaning differs**: `single_temp` values come from a single per-detector
current-pulse template, while `pulse_lib` values come from the per-pixel
pulse-shape library (see {ref}`hpge-psl-overview`). All fields are forwarded
from the corresponding `hit`-tier subtable and are `VectorOfVectors`
(variable-length per event).

| Field             | Units | Subtable                   | Description                                                                                                                                                                                             |
| ----------------- | ----- | -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `is_good`         |       | `psd`                      | Boolean. `True` if the PSD usability flag is valid in LEGEND-200 data.                                                                                                                                  |
| `is_valid_sim`    |       | `psd`                      | Boolean. `True` when the crystal metadata needed to model the detector is usable, so its simulated response can be trusted. Reshaped from the `hit`-tier `is_valid_sim` field (see {ref}`par-detinfo`). |
| `aoe`             |       | `single_temp`, `pulse_lib` | A/E classifier values forwarded from the `hit` tier.                                                                                                                                                    |
| `aoe_corr`        |       | `single_temp`, `pulse_lib` | Energy-corrected A/E values forwarded from the `hit` tier.                                                                                                                                              |
| `drift_time_amax` | ns    | `single_temp`, `pulse_lib` | Drift time at the maximum-current position, forwarded from the `hit` tier.                                                                                                                              |
| `has_aoe`         |       | `single_temp`, `pulse_lib` | Boolean. `True` if the A/E value is not `NaN` (i.e. PSD was computed).                                                                                                                                  |
| `is_single_site`  |       | `single_temp`, `pulse_lib` | Boolean single-site PSD flag forwarded from the `hit` tier.                                                                                                                                             |
| `is_bb_like`      |       | `pulse_lib`                | Boolean flag for $0\nu\beta\beta$-like single-site events, forwarded from `hit`.                                                                                                                        |
| `is_high_aoe`     |       | `pulse_lib`                | Boolean high-A/E flag forwarded from the `hit` tier.                                                                                                                                                    |

#### `geds/surface/` — energy depositions near the n+ surface

Present only with `store_surface_steps: true`. The fields of the `surface`
subtable of the [`hit` tier](hit-tier-output), for the same hits as
`geds/energy`.

| Field               | Type                               | Units | Description                                                                                |
| ------------------- | ---------------------------------- | ----- | ------------------------------------------------------------------------------------------ |
| `edep`              | `VectorOfVectors` (event/hit/step) | keV   | Energy of each step closer than `surface_depth_in_mm` to the n+ surface, one list per hit. |
| `distance_to_nplus` | `VectorOfVectors` (event/hit/step) | mm    | Distance of each of these steps to the n+ surface, in the same order as `edep`.            |
| `energy_bulk`       | `VectorOfVectors` (event/hit)      | keV   | Summed energy of all other steps of the hit.                                               |

### `spms/` — SiPM (LAr scintillation) array

Per-event arrays collecting SiPM data. All non-OFF channels are always present
in ascending UID order, even for events with no energy deposition in liquid
argon.

| Field          | Type              | Units | Description                                                                                                                                                                                                                                                                                     |
| -------------- | ----------------- | ----- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `rawid`        | `VectorOfVectors` | —     | SiPM channel UIDs in the simulated geometry. Named `rawid` to match LEGEND-200 data; equal to the data rawids only if the geometry was built from the same channel map. Always the full list of non-OFF channels per event.                                                                     |
| `energy`       | `VectorOfVectors` | —     | PE amplitudes per channel per event, filtered by the PE energy threshold. Nested variable-length array.                                                                                                                                                                                         |
| `time`         | `VectorOfVectors` | ns    | PE times per channel per event, relative to `trigger/timestamp`. Nested variable-length array matching `energy`. Unlike `spms/t0` in LEGEND-200 data, which counts from the start of the waveform.                                                                                              |
| `is_saturated` | `VectorOfVectors` | —     | Boolean SiPM saturation flag per channel. `True` if PE count exceeds threshold.                                                                                                                                                                                                                 |
| `hit_idx`      | `VectorOfVectors` | —     | Row index in the `opt`-tier table for lookback (the first row, if the channel has several in the event). Set to `-1` for channels without an `opt` row in the event.                                                                                                                            |
| `expected_pes` | `VectorOfVectors` | —     | _(optional)_ Expected number of photoelectrons per channel at unit channel efficiency, before the PE-per-hit cap and the PE threshold. Same order as `rawid`; `0` for channels without an `opt` row in the event. Present only when `store_expected_pes` is enabled in the `opt` tier settings. |
| `energy_sum`   | `Array`           | —     | Total PE energy summed over all channels and all PEs. Includes the random-coincidence PEs if `add_random_coincidences` is enabled. Scalar per event.                                                                                                                                            |
| `multiplicity` | `Array`           | —     | Number of SiPM channels with at least one PE. Includes the random-coincidence PEs if `add_random_coincidences` is enabled. Scalar per event.                                                                                                                                                    |
| `rc_energy`    | `VectorOfVectors` | —     | _(optional)_ Random-coincidence PE amplitudes taken from LEGEND-200 data (see `random_coincidence_mode`), in the `rawid` channel order. Channels are matched by name through the rawids of the source run. Present only when `add_random_coincidences` is enabled.                              |
| `rc_time`      | `VectorOfVectors` | ns    | _(optional)_ Random-coincidence PE times, on the same time axis as `time`. Present only when `add_random_coincidences` is enabled.                                                                                                                                                              |

### `coincident/` — detector coincidence flags

| Field  | Type    | Units | Description                                                                          |
| ------ | ------- | ----- | ------------------------------------------------------------------------------------ |
| `geds` | `Array` | —     | Boolean. `True` if the HPGe multiplicity is greater than zero.                       |
| `spms` | `Array` | —     | Boolean LAr veto flag. `True` if `spms/multiplicity >= 4` or `spms/energy_sum >= 4`. |

## Time-coincidence map (TCM)

Every `hit`, `opt`, and `evt` tier file contains a `/tcm` table
(time-coincidence map) that maps physics events to the individual detector-level
table rows that belong to them. The TCM is built by grouping hits that share the
same `evtid` and whose `t0` values fall within a 10 µs coincidence window
(matching the _remage_ built-in TCM settings).

The TCM table has two fields, both `VectorOfVectors` (one inner list per event):

| Field          | Type              | Description                                                                                                  |
| -------------- | ----------------- | ------------------------------------------------------------------------------------------------------------ |
| `table_key`    | `VectorOfVectors` | Detector UID for each hit in the event. Identifies which `/hit/{detector}/` table the hit belongs to.        |
| `row_in_table` | `VectorOfVectors` | Row index into the corresponding detector table. Together with `table_key`, uniquely locates each hit entry. |

In the `hit` and `opt` tiers, the TCM indexes into the detector tables within
the same file. In the `evt` tier, the TCM is a _unified_ version that merges the
`hit` and `opt` TCMs, so that a single TCM entry references hits across both
HPGe and SiPM detector tables.

(pdf-tier)=

## `pdf` tier — probability density functions

The `pdf` tier reads the event-level data from the `cvt` tier and bins it into
energy histograms. These histograms represent the probability density functions
(PDFs) used as inputs to spectral fitting analyses. The output is a single LH5
file containing a set of histograms, each corresponding to a different event
selection and detector group, and a scalar recording the total number of
simulated primary events.

The histograms split the events into independent categories by multiplicity, LAr
veto and pulse-shape discrimination (PSD). All histograms carry a `description`
attribute in the LH5 attrs.

### Root-level fields

| Field           | Type     | Units | Description                                                                                                                                                                                                                             |
| --------------- | -------- | ----- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `nr_sim_events` | `Scalar` | —     | Total number of simulated primary events, read from the `number_of_simulated_events` scalar that _remage_ stores in each `stp` file and that is summed over all jobs at the `cvt` tier. Used to normalise PDFs to physical event rates. |

### `pdf/` histogram struct

The 1D histograms are stored at paths of the form

```text
pdf/mul1/<lar_category>/<aoe_category>[/<fccd_category>]/<group>
```

Only multiplicity-1 events are histogrammed in 1D: exactly one detector fired
(`geds.multiplicity == 1`), and it is ON (not AC or OFF, i.e.
`geds.quality.is_good_channel` is `True`). Hits below the `evt` tier energy
threshold `geds_energy_thr_kev` are ignored. The category levels below `mul1`
split these events by the LAr veto, by A/E and, optionally, by the energy
deposited in the FCCD layer. Every event falls in exactly one category of each
level, so the histograms do not overlap.

The last path element is the detector group. Each group is a `Histogram` of the
energies (in keV, 1 keV bins from 0 to 6000 keV) of the hits in the detectors of
that group. The categories select whole events, independent of the group. The
groups are configured with the `detector_groups` setting (see
{ref}`pdf-tier-settings`), and the `all` group with every detector is always
present. For example, with `detector_groups: {icpc: "V.*", bege: "B.*"}`, each
category contains `icpc`, `bege` and `all`.

Example with the pulse-shape library and the `all` group only:

```text
pdf/
├── mul2/all
└── mul1/
    ├── is_spms_coincident/
    │   ├── is_high_aoe/all
    │   ├── is_ss_aoe/all
    │   ├── is_low_aoe/all
    │   └── is_invalid_aoe/all
    └── not_spms_coincident/
        ├── is_high_aoe/all
        ├── is_ss_aoe/all
        ├── is_low_aoe/all
        └── is_invalid_aoe/all
```

#### Categories

The categories are computed by the functions in {mod}`legendsimflow.pdf`, which
other scripts can use too.

| Level | Category              | Selection                                                                                     |
| ----- | --------------------- | --------------------------------------------------------------------------------------------- |
| LAr   | `is_spms_coincident`  | light seen by the SiPMs (`coincident.spms` is `True`), i.e. the event fails the LAr veto      |
| LAr   | `not_spms_coincident` | no light seen by the SiPMs (`coincident.spms` is `False`), i.e. the event passes the LAr veto |
| A/E   | `is_high_aoe`         | A/E above the high-side cut (`is_high_aoe`). Pulse-shape library only                         |
| A/E   | `is_ss_aoe`           | single-site: A/E above the low-side cut (`is_single_site`) and below the high-side cut        |
| A/E   | `is_low_aoe`          | multi-site: A/E below the low-side cut                                                        |
| A/E   | `is_invalid_aoe`      | no A/E available: `psd.is_good` or `has_aoe` is `False`                                       |
| FCCD  | `is_fccd_tagged`      | more than `fccd_tag_energy_thr_kev` deposited in the FCCD layer                               |
| FCCD  | `not_fccd_tagged`     | at most `fccd_tag_energy_thr_kev` deposited in the FCCD layer                                 |

The A/E flags are taken from the pulse-shape library model (`psd.pulse_lib`)
when the `hit` tier runs with `simulate_psd_with_psl`, and otherwise from the
single-template model (`psd.single_temp`) when it runs with `simulate_psd` (see
{ref}`hit-tier-settings`). The single-template model has no high-side cut, so
`is_high_aoe` is absent. This gives 8 categories with the pulse-shape library
and 6 with the single template, twice as many with the FCCD level.

`is_invalid_aoe` holds events with a detector that is ON but whose PSD response
could not be simulated (e.g. because it is not included in the simulation
model), or that has no valid PSD in data. An analysis applying an A/E cut treats
these events as background.

The FCCD level is present only when `fccd_tag_energy_thr_kev` is set (see
{ref}`pdf-tier-settings`). The energy in the FCCD layer is the sum of the
depositions in `geds/surface/edep` closer to the n+ surface than the FCCD of the
detector, as deposited, before the dead-layer model. The FCCD is read from the
detector metadata. It is limited to `surface_depth_in_mm`, the depth of the
stored depositions, with a warning.

A level is omitted when the tier settings say its data was not produced: the LAr
level with `skip_opt` in the `evt` tier, the A/E level when the `hit` tier runs
with neither `simulate_psd` nor `simulate_psd_with_psl`. For example, with
`skip_opt` the paths are `pdf/mul1/<aoe_category>/<group>`.

:::{tip}

Any other selection is a sum of categories, which
{func}`legendsimflow.pdf.read_pdf` reads directly. A path selects every category
below it, and `*` matches any category at its level:

```python
from legendsimflow.pdf import read_pdf

read_pdf("mul1", "pdf.lh5")  # all multiplicity-1 events
read_pdf("mul1/not_spms_coincident", "pdf.lh5")  # after the LAr veto
read_pdf("mul1/*/is_ss_aoe", "pdf.lh5")  # after the A/E cuts only
read_pdf("mul1/not_spms_coincident/is_ss_aoe", "pdf.lh5", group="icpc")
```

:::

#### 2D histograms

| Key        | Description                                                                                                                                                                   |
| ---------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mul2/all` | Multiplicity-2 events with exactly two ON detectors fired. A 2D histogram with axes (E_low, E_high), where E_low ≤ E_high are the two hit energies sorted in ascending order. |

`mul2` is **not** split by detector group: only the `all` group is present.
Per-channel or per-pair 2-D PDFs are out of scope for the current
implementation.

### Spectra on a grid of parameters

{func}`legendsimflow.pdf.build_pdf_grid` evaluates a PDF, a function of the
parameter values that returns a histogram, on an evenly spaced grid of
parameters. [Refit.jl](https://github.com/legend-exp/Refit.jl) reads the result
as a generalised histogram and interpolates between the grid points. The
parameters can be anything the PDF depends on; the dead-layer model below is one
example.

With `store_surface_steps: true`, the PDF can depend on the dead-layer model.
For the hits of one detector (`in_det`) in the `evt/geds` table of a `cvt` file
(`geds`), with nominal FCCD `fccd` (mm) and dead layer fraction `dlf`:

```python
surface = ak.to_packed(geds.surface[in_det])  # copy the hits of the detector
energy = ak.flatten(geds.energy[in_det])
dist = ak.flatten(surface.distance_to_nplus)
edep = ak.flatten(surface.edep)
act_nominal = piecewise_linear_activeness(dist, fccd, dlf)


def pdf(pars):
    act = piecewise_linear_activeness(dist, pars["fccd"], pars["dlf"])
    # nominal (smeared) energy plus the change of the collected energy
    e = energy + ak.sum(edep * (act - act_nominal), axis=-1)
    return hist.new.Reg(6000, 0, 6000).Double().fill(e[e > 25])  # keV threshold


grid = build_pdf_grid(pdf, {"fccd": (0.5, 1.5, 0.1), "dlf": (0, 0.9, 0.3)})
```

Cuts that depend on the energy, such as the threshold or the multiplicity, are
applied to the recomputed energies inside `pdf`.

:::{warning}

The `evt` tier drops the hits below the energy threshold (`geds_energy_thr_kev`,
see {ref}`evt-tier-settings-meta`), computed with the nominal dead-layer model,
also with `store_surface_steps: true`. A hit below the threshold for the nominal
model can pass it for a thinner FCCD, but it is missing from the `cvt` file. The
recomputed PDF is therefore biased near the threshold, and so is the
multiplicity of the events.

:::
