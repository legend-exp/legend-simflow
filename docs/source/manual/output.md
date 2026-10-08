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
modeling flags for **every deployed germanium** detector:

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

(opt-tier)=

## `opt` tier — optical (SiPM) post-processing

The `opt` tier is at the same conceptual level as the `hit` tier: it performs
detector-wise post-processing, but for SiPMs instead of HPGe detectors. It
applies the optical map convolution and photoelectron (PE) response models to
the scintillator output from the `stp` tier. When using per-SiPM optical maps, a
separate table is written for each SiPM channel under `/hit/{sipm_name}/`; when
using a single summed map (the default), all SiPM channels are aggregated into a
single `/hit/spms/` table. Each row corresponds to an `stp`-tier hit entry
(identified by `evtid`) — not to a single physics event: a liquid argon hit with
the optical map, a hit of the SiPM itself with `light_source: tracked_photons`
(see {ref}`opt-tier-settings`). In the latter case every SiPM of the geometry
has a table, empty if no photon reached it.

### Inherited fields

| Field   | Type    | Units | Description                                                                                                                                        |
| ------- | ------- | ----- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| `evtid` | `Array` | —     | Event identifier, shared across all detectors hit in the same event.                                                                               |
| `t0`    | `Array` | ns    | Time of the first energy deposition in the LAr hit (first photon in the SiPM hit with tracked photons), relative to the start of the Geant4 event. |

### Added fields

| Field          | Type              | Units | Description                                                                                                                                                                            |
| -------------- | ----------------- | ----- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `dt`           | `VectorOfVectors` | ns    | Photoelectron arrival times relative to the hit `t0`, after resolution smearing and photoelectron clustering to simulate the timing resolution of the SiPM. Variable-length per row.   |
| `energy`       | `VectorOfVectors` | —     | Photoelectron amplitudes (relative units), after PE resolution smearing. Variable-length array matching `dt`.                                                                          |
| `is_saturated` | `Array`           | —     | Boolean flag. `True` when the number of detected photoelectrons exceeds a maximum PE-per-hit cap, indicating SiPM saturation.                                                          |
| `expected_pes` | `Array`           | —     | _(optional)_ Expected number of photoelectrons per row at unit channel efficiency, before the PE-per-hit cap. Present only when `store_expected_pes` is enabled, with the optical map. |
| `period`       | `Array`           | —     | Data-taking period number extracted from the run identifier (numeric encoding).                                                                                                        |
| `run`          | `Array`           | —     | Data-taking run number extracted from the run identifier (numeric encoding).                                                                                                           |
| `usability`    | `Array`           | —     | Encoded SiPM channel usability status for this run. Decode with {func}`legendsimflow.metadata.decode_usability`.                                                                       |

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

| Field       | Type    | Units | Description                                                                                                                                                                                                                                                                                       |
| ----------- | ------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `evtid`     | `Array` | —     | Event identifier.                                                                                                                                                                                                                                                                                 |
| `period`    | `Array` | —     | Data-taking period number.                                                                                                                                                                                                                                                                        |
| `run`       | `Array` | —     | Data-taking run number.                                                                                                                                                                                                                                                                           |
| `timestamp` | `Array` | ns    | `t0` of the first HPGe hit, or of the first SiPM hit when the `hit` tier is skipped (`skip_hit: true`), relative to the start of the simulated Geant4 event. Used as the event's timestamp. Stored in double precision, because decays late in a chain happen hours after the start of the event. |

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

### `spms/` — SiPM (LAr scintillation) array

Per-event arrays collecting SiPM data. All non-OFF channels are always present
in ascending UID order, even for events with no energy deposition in liquid
argon.

| Field          | Type              | Units | Description                                                                                                                                                                                                                                                                                     |
| -------------- | ----------------- | ----- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `rawid`        | `VectorOfVectors` | —     | SiPM channel UIDs, matching the channel identifiers used in LEGEND-200 data. Always the full list of non-OFF channels per event.                                                                                                                                                                |
| `energy`       | `VectorOfVectors` | —     | PE amplitudes per channel per event, filtered by the PE energy threshold. Nested variable-length array.                                                                                                                                                                                         |
| `time`         | `VectorOfVectors` | ns    | PE times per channel per event, relative to `trigger/timestamp`. Nested variable-length array matching `energy`. Unlike `spms/t0` in LEGEND-200 data, which counts from the start of the waveform.                                                                                              |
| `is_saturated` | `VectorOfVectors` | —     | Boolean SiPM saturation flag per channel. `True` if PE count exceeds threshold.                                                                                                                                                                                                                 |
| `hit_idx`      | `VectorOfVectors` | —     | Row index in the `opt`-tier table for lookback (the first row, if the channel has several in the event). Set to `-1` for channels without an `opt` row in the event.                                                                                                                            |
| `expected_pes` | `VectorOfVectors` | —     | _(optional)_ Expected number of photoelectrons per channel at unit channel efficiency, before the PE-per-hit cap and the PE threshold. Same order as `rawid`; `0` for channels without an `opt` row in the event. Present only when `store_expected_pes` is enabled in the `opt` tier settings. |
| `energy_sum`   | `Array`           | —     | Total PE energy summed over all channels and all PEs. Scalar per event.                                                                                                                                                                                                                         |
| `multiplicity` | `Array`           | —     | Number of SiPM channels with at least one detected PE. Scalar per event.                                                                                                                                                                                                                        |
| `rc_energy`    | `VectorOfVectors` | —     | _(optional)_ Random-coincidence PE amplitudes from forced-trigger data. Present only when `add_random_coincidences` is enabled.                                                                                                                                                                 |
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

The histograms apply a sequence of analysis cuts — multiplicity, LAr
anti-coincidence, and pulse-shape discrimination — to produce PDFs for the most
common LEGEND-200 analysis channels. All histograms carry a `description`
attribute in the LH5 attrs.

### Root-level fields

| Field           | Type     | Units | Description                                                                                                                                                                                                                             |
| --------------- | -------- | ----- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `nr_sim_events` | `Scalar` | —     | Total number of simulated primary events, read from the `number_of_simulated_events` scalar that _remage_ stores in each `stp` file and that is summed over all jobs at the `cvt` tier. Used to normalise PDFs to physical event rates. |

### `pdf/` histogram struct

The 1D histograms are stored at paths of the form

```text
pdf/<level_1>/<level_2>/.../<level_n>/<group>
```

Each level is one event selection, and a histogram contains the events that
satisfy all levels on its path. A level is either `<cut>`, for events passing
the cut, or `not_<cut>`, for events failing it. The levels always appear in the
order `mul1`, `lar`, PSD cut. A cut absent from a path is not applied.

The last path element is the detector group. Each group is a `Histogram` of the
energies (in keV, 1 keV bins from 0 to 6000 keV) of the hits in the detectors of
that group. The cuts select whole events, independent of the group. The groups
are configured with the `detector_groups` setting (see
{ref}`pdf-tier-settings`), and the `all` group with every detector is always
present. For example, with `detector_groups: {icpc: "V.*", bege: "B.*"}`, the
output contains `pdf/mul1/lar/icpc`, `pdf/mul1/lar/bege` and `pdf/mul1/lar/all`.
A level holds the group histograms of its selection next to the deeper levels,
so group names cannot be cut levels (`lar`, `not_lar`, `aoe_st`, ...).

Example, with the `all` group only:

```text
pdf/
├── hit/all
├── mul2
└── mul1/
    ├── all
    ├── lar/
    │   ├── all
    │   └── aoe_psl/all
    ├── not_lar/all
    ├── aoe_psl/all
    └── not_aoe_psl/all
```

`pdf/mul1/lar/aoe_psl/all` contains multiplicity-1 events with no light in the
LAr instrumentation that pass both A/E cuts.

#### Levels

All selections require every hit in the event to be in an ON detector (not AC or
OFF), i.e. `geds.quality.is_good_channel` is `True`.

| Level  | Selection                                                                                                                                          |
| ------ | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| `hit`  | All HPGe energy deposits, no multiplicity requirement. Top level only, no further levels below.                                                    |
| `mul1` | Multiplicity-1 events: exactly one ON detector fired (`geds.multiplicity == 1`). All other levels are below `mul1`.                                |
| `lar`  | LAr anti-coincidence: no light seen by the SiPMs (`coincident.spms` is `False`). `not_lar` selects events with light. Present only with SiPM data. |
| PSD    | PSD cuts, see the table below. Each requires `psd.is_good` and the `has_aoe` flag of its A/E model for all hits, plus the cut condition.           |

| PSD level      | A/E model           | Condition on all hits                                    |
| -------------- | ------------------- | -------------------------------------------------------- |
| `aoe_st`       | single template     | `psd.single_temp.is_single_site` (low-side A/E cut)      |
| `aoe_psl_low`  | pulse-shape library | `psd.pulse_lib.is_single_site` (low-side A/E cut)        |
| `aoe_psl_high` | pulse-shape library | not `psd.pulse_lib.is_high_aoe` (high-side A/E cut)      |
| `aoe_psl`      | pulse-shape library | `psd.pulse_lib.is_bb_like` (low- and high-side A/E cuts) |

The single-template levels are present when the `hit` tier runs with
`simulate_psd`, the pulse-shape library levels when it runs with
`simulate_psd_with_psl` (see {ref}`hit-tier-settings`).

The following paths are written, for each PSD level `<aoe>` available:

| Path                           | Selection                                        |
| ------------------------------ | ------------------------------------------------ |
| `hit`                          | all HPGe energy deposits                         |
| `mul1`                         | multiplicity 1                                   |
| `mul1/lar`, `mul1/not_lar`     | multiplicity 1, passing / failing the LAr veto   |
| `mul1/<aoe>`, `mul1/not_<aoe>` | multiplicity 1, passing / failing the PSD cut    |
| `mul1/lar/<aoe>`               | multiplicity 1, passing the LAr veto and PSD cut |

:::{warning}

A PSD cut has three outcomes. Events with `psd.is_good` or `has_aoe` set to
`False` for any hit pass neither `<aoe>` nor `not_<aoe>`. This is the case when
a detector is ON with valid PSD in data but its PSD response could not be
simulated (e.g. because it is not included in the simulation model). These
events are treated as background: rather than keeping events we cannot
characterise, we cut them. As a consequence, `mul1/<aoe>` and `mul1/not_<aoe>`
do not add up to `mul1`.

:::

#### 2D histograms

| Key    | Description                                                                                                                                                                   |
| ------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mul2` | Multiplicity-2 events with exactly two ON detectors fired. A 2D histogram with axes (E_low, E_high), where E_low ≤ E_high are the two hit energies sorted in ascending order. |

`mul2` is a single global histogram and is **not** split by detector group.
Per-channel or per-pair 2-D PDFs are out of scope for the current
implementation.
