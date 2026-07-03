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
carry `usability` alone, while `psd_usability`, `crystal_metadata_usability`,
`is_modelable`, and `operational_voltage_in_V` are germanium-only.

The [`cache_detector_usabilities`](../api/snakemake_rules.md) rule writes the
data-quality flags for **all** deployed detectors:

| File                              | Applies to  | Value                | Description                                                                                                                                       |
| --------------------------------- | ----------- | -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| `usability.yaml`                  | geds + SiPM | `on`, `off`, `ac`, … | Analysis usability from the channel-map status (`analysis.usability`).                                                                            |
| `psd_usability.yaml`              | geds        | `valid`, …           | PSD usability from `analysis.psd.status.low_aoe`; defaults to `valid` when the field is absent.                                                   |
| `crystal_metadata_usability.yaml` | geds        | `valid`, …, `null`   | Usability of the crystal metadata required for modeling (the crystal-slice `status`). `null` when the information is unavailable in the metadata. |

The [`cache_modelable_hpges`](../api/snakemake_rules.md) checkpoint writes the
modeling flags for **every deployed germanium** detector:

| File                            | Value            | Description                                                                                                                                     |
| ------------------------------- | ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `is_modelable.yaml`             | `true` / `false` | Whether the detector is suitable for drift-time-map and current-pulse modeling. See {ref}`hpge-modeling-criteria` for the eligibility criteria. |
| `operational_voltage_in_V.yaml` | integer / `null` | Operational bias voltage read from the parameters database. `null` for detectors with no bias (e.g. off detectors).                             |

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

The `hit` tier can add HPGe pulse-shape-discrimination (PSD) fields in two
optional subtables, selected by the metadata settings in
{ref}`hit-tier-settings`:

- `psd` — single-template A/E simulation, present only when `simulate_psd: True`
  (the default).
- `psd_psl` — pulse-shape-library (PSL) based simulation, present only when
  `simulate_psd_with_psl: True` (see {ref}`hpge-psl-overview`).

:::{note}

The two subtables hold the same-named A/E fields, but their **meaning differs**:
the `psd` values come from a single per-detector current-pulse template, while
the `psd_psl` values come from the per-pixel pulse-shape library (see
{ref}`hpge-psl-overview`).

:::

| Field             | Type    | Units | Subtable         | Description                                                                                                                                                                                                                                                                                            |
| ----------------- | ------- | ----- | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `drift_time_amax` | `Array` | ns    | `psd`, `psd_psl` | Drift time at the maximum-current (A) position of the simulated current pulse. Set to `NaN` when no drift-time map or current-pulse model is available for the detector.                                                                                                                               |
| `aoe_raw`         | `Array` |       | `psd`, `psd_psl` | Raw A/E value: maximum current amplitude divided by energy. The maximum current (A) is obtained from the simulated current pulse and includes electronic noise effects.                                                                                                                                |
| `aoe_corr`        | `Array` |       | `psd`, `psd_psl` | Energy-corrected A/E value, obtained by correcting `aoe_raw` for the observed energy dependence.                                                                                                                                                                                                       |
| `aoe`             | `Array` |       | `psd`, `psd_psl` | A/E classifier value: `(aoe_corr - 1) / aoe_resolution`. Used for pulse-shape discrimination (PSD). Set to `NaN` when PSD simulation is not available (i.e., when the drift-time map or current-pulse model is missing). This is distinct from the usability flags that track LEGEND-200 data quality. |
| `is_single_site`  | `Array` |       | `psd`, `psd_psl` | Boolean PSD flag. `True` when the A/E classifier `aoe` exceeds the lower single-site cut `psdcuts.aoe.low_side` (extracted from LEGEND-200 data).                                                                                                                                                      |
| `is_bb_like`      | `Array` |       | `psd_psl`        | Boolean PSD flag. `True` for $0\nu\beta\beta$-like single-site events: `aoe` above `psdcuts.aoe.low_side` and not above the upper cut `psdcuts.aoe.high_side`.                                                                                                                                         |
| `is_high_aoe`     | `Array` |       | `psd_psl`        | Boolean PSD flag. `True` when `aoe` exceeds the upper cut `psdcuts.aoe.high_side` (high-A/E events, e.g. surface or $\alpha$).                                                                                                                                                                         |

## `opt` tier — optical (SiPM) post-processing

The `opt` tier is at the same conceptual level as the `hit` tier: it performs
detector-wise post-processing, but for SiPMs instead of HPGe detectors. It
applies the optical map convolution and photoelectron (PE) response models to
the scintillator output from the `stp` tier. When using per-SiPM optical maps, a
separate table is written for each SiPM channel under `/hit/{sipm_name}/`; when
using a single summed map (the default), all SiPM channels are aggregated into a
single `/hit/spms/` table. Each row corresponds to an `stp`-tier hit entry
(identified by `evtid`) — not to a single physics event.

### Inherited fields

| Field   | Type    | Units | Description                                                          |
| ------- | ------- | ----- | -------------------------------------------------------------------- |
| `evtid` | `Array` | —     | Event identifier, shared across all detectors hit in the same event. |

### Added fields

| Field          | Type              | Units | Description                                                                                                                                             |
| -------------- | ----------------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `time`         | `VectorOfVectors` | ns    | Photoelectron hit times, after resolution smearing and photoelectron clustering to simulate the timing resolution of the SiPM. Variable-length per row. |
| `energy`       | `VectorOfVectors` | —     | Photoelectron amplitudes (relative units), after PE resolution smearing. Variable-length array matching `time`.                                         |
| `is_saturated` | `Array`           | —     | Boolean flag. `True` when the number of detected photoelectrons exceeds a maximum PE-per-hit cap, indicating SiPM saturation.                           |
| `period`       | `Array`           | —     | Data-taking period number extracted from the run identifier (numeric encoding).                                                                         |
| `run`          | `Array`           | —     | Data-taking run number extracted from the run identifier (numeric encoding).                                                                            |
| `usability`    | `Array`           | —     | Encoded SiPM channel usability status for this run. Decode with {func}`legendsimflow.metadata.decode_usability`.                                        |

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

| Field       | Type    | Units | Description                                                                                                                 |
| ----------- | ------- | ----- | --------------------------------------------------------------------------------------------------------------------------- |
| `evtid`     | `Array` | —     | Event identifier.                                                                                                           |
| `period`    | `Array` | —     | Data-taking period number.                                                                                                  |
| `run`       | `Array` | —     | Data-taking run number.                                                                                                     |
| `timestamp` | `Array` | ns    | Delta time of the first HPGe hit (`t0`) relative to the start of the simulated Geant4 event. Used as the event's timestamp. |

### `geds/` — HPGe detector array

Per-event arrays collecting HPGe hits that pass the energy threshold (25 keV)
and are from non-OFF detectors.

| Field             | Type              | Units | Description                                                                                                             |
| ----------------- | ----------------- | ----- | ----------------------------------------------------------------------------------------------------------------------- |
| `energy`          | `VectorOfVectors` | keV   | Hit energies from ON and AC detectors above threshold. Variable-length per event.                                       |
| `energy_sum`      | `Array`           | keV   | Summed energy from ON detectors only (excludes AC). Scalar per event.                                                   |
| `rawid`           | `VectorOfVectors` | —     | Detector channel UID for each hit, matching the channel identifiers used in LEGEND-200 data. Variable-length per event. |
| `hit_idx`         | `VectorOfVectors` | —     | Row index in the `hit`-tier table, for looking up additional hit-level fields. Variable-length per event.               |
| `is_good_channel` | `VectorOfVectors` | —     | Boolean. `True` if the detector usability is ON (not AC or OFF). Variable-length per event.                             |
| `multiplicity`    | `Array`           | —     | Number of HPGe hits above threshold per event. Scalar per event.                                                        |

#### `geds/psd/` and `geds/psd_psl` — PSD fields

These two subtables are optional and gated by the `hit`-tier settings
({ref}`hit-tier-settings`): `geds/psd` is present when `simulate_psd: True` and
`geds/psd_psl` when `simulate_psd_with_psl: True`. Each carries the same-named
fields forwarded from the corresponding `hit`-tier subtable (`psd` is
single-template, `psd_psl` is PSL based; see {ref}`hpge-psl-overview`), so a
field's meaning depends on which subtable it is in. All fields are
`VectorOfVectors` (variable-length per event).

| Field             | Units | Subtable         | Description                                                                                                                                                                                             |
| ----------------- | ----- | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `is_good`         |       | `psd`            | Boolean. `True` if the PSD usability flag is valid in LEGEND-200 data.                                                                                                                                  |
| `is_valid_sim`    |       | `psd`            | Boolean. `True` when the crystal metadata needed to model the detector is usable, so its simulated response can be trusted. Reshaped from the `hit`-tier `is_valid_sim` field (see {ref}`par-detinfo`). |
| `aoe`             |       | `psd`, `psd_psl` | A/E classifier values forwarded from the `hit` tier.                                                                                                                                                    |
| `aoe_corr`        |       | `psd`, `psd_psl` | Energy-corrected A/E values forwarded from the `hit` tier.                                                                                                                                              |
| `drift_time_amax` | ns    | `psd`, `psd_psl` | Drift time at the maximum-current position, forwarded from the `hit` tier.                                                                                                                              |
| `has_aoe`         |       | `psd`, `psd_psl` | Boolean. `True` if the A/E value is not `NaN` (i.e. PSD was computed).                                                                                                                                  |
| `is_single_site`  |       | `psd`, `psd_psl` | Boolean single-site PSD flag forwarded from the `hit` tier.                                                                                                                                             |
| `is_bb_like`      |       | `psd_psl`        | Boolean flag for $0\nu\beta\beta$-like single-site events, forwarded from `hit`.                                                                                                                        |
| `is_high_aoe`     |       | `psd_psl`        | Boolean high-A/E flag forwarded from the `hit` tier.                                                                                                                                                    |

### `spms/` — SiPM (LAr scintillation) array

Per-event arrays collecting SiPM data. All non-OFF channels are always present
in ascending UID order, even for events with no energy deposition in liquid
argon.

| Field          | Type              | Units | Description                                                                                                                      |
| -------------- | ----------------- | ----- | -------------------------------------------------------------------------------------------------------------------------------- |
| `rawid`        | `VectorOfVectors` | —     | SiPM channel UIDs, matching the channel identifiers used in LEGEND-200 data. Always the full list of non-OFF channels per event. |
| `energy`       | `VectorOfVectors` | —     | PE amplitudes per channel per event, filtered by the PE energy threshold. Nested variable-length array.                          |
| `time`         | `VectorOfVectors` | ns    | PE hit times per channel per event. Nested variable-length array matching `energy`.                                              |
| `is_saturated` | `VectorOfVectors` | —     | Boolean SiPM saturation flag per channel. `True` if PE count exceeds threshold.                                                  |
| `hit_idx`      | `VectorOfVectors` | —     | Row index in the `opt`-tier table for lookback. Set to `-1` for events with no LAr energy deposition.                            |
| `energy_sum`   | `Array`           | —     | Total PE energy summed over all channels and all PEs. Scalar per event.                                                          |
| `multiplicity` | `Array`           | —     | Number of SiPM channels with at least one detected PE. Scalar per event.                                                         |
| `rc_energy`    | `VectorOfVectors` | —     | _(optional)_ Random-coincidence PE amplitudes from forced-trigger data. Present only when `add_random_coincidences` is enabled.  |
| `rc_time`      | `VectorOfVectors` | ns    | _(optional)_ Random-coincidence PE times. Present only when `add_random_coincidences` is enabled.                                |

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

### `pdf/` — histogram struct

The 1-D histograms are organised cut-first, then by detector group:
`pdf/<cut>/<group>`. Each leaf is an lgdo `Histogram`. The detector groups are
configured via the `detector_groups` setting (see {ref}`pdf-tier-settings`); the
`all` group containing every detector is always present.

For example, with `detector_groups: {icpc: "V.*", bege: "B.*"}`, the output
contains `pdf/hit/icpc`, `pdf/hit/bege`, and `pdf/hit/all`, and similarly for
every other cut.

#### 1D histograms

Each of the following cuts produces one `Histogram` per detector group. The
`good_channel_mask` applied before all cuts requires every channel in the event
to be an ON detector (not AC or OFF). Per-group filtering restricts which
detector energies are accumulated into the histogram; the event-level cuts
themselves are unchanged and applied globally.

| Cut           | Description                                                                                                                                                                                                                    |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `hit`         | All individual HPGe energy deposits in ON-channel events, with no multiplicity requirement.                                                                                                                                    |
| `mul`         | Multiplicity-1 events: exactly one ON detector fired (`geds.multiplicity == 1`).                                                                                                                                               |
| `mul_lar`     | Multiplicity-1 events passing the LAr anti-coincidence cut. Events are vetoed when `coincident.spms` is `True` (SiPMs detected scintillation light in liquid argon). Present only when SiPM data is available.                 |
| `mul_psd`     | Multiplicity-1 events passing the PSD single-site cut. Requires `psd.is_good`, `psd.has_aoe`, and `psd.is_single_site` for all hits. Events where PSD is not valid or not simulated are classified as background and excluded. |
| `mul_lar_psd` | Multiplicity-1 events passing both the LAr anti-coincidence and PSD single-site cuts (combination of `mul_lar` and `mul_psd`). Present only when SiPM data is available.                                                       |

:::{warning}

When a detector is ON with valid PSD in the data but its PSD response could not
be simulated (e.g. because it is not included in the simulation model), the
corresponding events will have `psd.has_aoe = False` in the `cvt` tier. Such
events are treated as background and excluded from `mul_psd` and `mul_lar_psd`.
This is a conservative choice: rather than keeping events we cannot
characterise, we cut them. The `fail/psd` histogram does **not** include these
events either, since it is restricted to events where both `psd.is_good = True`
and `psd.has_aoe = True`.

:::

#### 2D histograms

| Key    | Description                                                                                                                                                                   |
| ------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mul2` | Multiplicity-2 events with exactly two ON detectors fired. A 2D histogram with axes (E_low, E_high), where E_low ≤ E_high are the two hit energies sorted in ascending order. |

`mul2` is a single global histogram and is **not** split by detector group.
Per-channel or per-pair 2-D PDFs are out of scope for the current
implementation.

### `pdf/fail/` — cut-failure histograms

The `fail/` sub-struct contains histograms for multiplicity-1 events that are
explicitly rejected by a cut, providing a way to characterise the vetoed
background. Like the pass histograms, each cut contains one `Histogram` per
detector group (`pdf/fail/<cut>/<group>`).

| Cut   | Description                                                                                                                                                                 |
| ----- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `lar` | Multiplicity-1 events failing the LAr veto (`coincident.spms == True`). Present only when SiPM data is available.                                                           |
| `psd` | Multiplicity-1 events with valid PSD (`psd.is_good == True`) that fail the single-site cut (`psd.is_single_site == False`). Events without valid PSD are not included here. |
