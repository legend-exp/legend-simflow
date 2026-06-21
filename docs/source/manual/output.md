# Output tier fields

This page documents the fields (columns) produced by each post-processing tier
of the Simflow. All output files use the
[LEGEND HDF5 (LH5) format](https://legend-exp.github.io/legend-data-format-specs/dev/hdf5/).

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

| Field           | Type    | Units | Description                                                                                                                                                                                                 |
| --------------- | ------- | ----- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `energy`        | `Array` | keV   | Reconstructed energy after smearing with the detector energy resolution. Computed from the sum of active energy depositions (weighted by the dead-layer activeness model).                                  |
| `period`        | `Array` | —     | Data-taking period number extracted from the run identifier (numeric encoding).                                                                                                                             |
| `run`           | `Array` | —     | Data-taking run number extracted from the run identifier (numeric encoding).                                                                                                                                |
| `usability`     | `Array` | —     | Encoded detector usability status for this run (e.g. `on`, `off`, `ac`). Decode with {func}`legendsimflow.metadata.decode_usability`. See the detector status flags in `legend-metadata/datasets/statuses`. |
| `psd_usability` | `Array` | —     | Encoded PSD usability flag (e.g. `valid`). Indicates whether PSD parameters are valid in LEGEND-200 data for this detector and run. Decode with {func}`legendsimflow.metadata.decode_psd_usability`.        |

In addition, several PSD based fields can be added, these are in either the
subtable `psd`, for the single template based A/E simulation (present if
`simulate_psd` is `True` in the hit tier setting file), or `psd_psl` for the
pulse shape library based one (present if `simulate_psd_with_psl` is `True` in
the hit tier settings file). Each subtable contains the following fields:

| Field             | Type    | Units | Description                                                                                                                                                                                                                                                                                            |
| ----------------- | ------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `drift_time_amax` | `Array` | ns    | Drift time at the maximum-current (A) position of the simulated current pulse. Set to `NaN` when no drift-time map or current-pulse model is available for the detector.                                                                                                                               |
| `aoe_raw`         | `Array` | —     | Raw A/E value: maximum current amplitude divided by energy. The maximum current (A) is obtained from the simulated current pulse, constructed from individual hit drift times and a current-pulse model, and includes electronic noise effects.                                                        |
| `aoe_corr`        | `Array` | —     | Energy-corrected A/E value, obtained by correcting `aoe_raw` for the observed energy dependence.                                                                                                                                                                                                       |
| `aoe`             | `Array` | —     | A/E classifier value: `(aoe_corr - 1) / aoe_resolution`. Used for pulse-shape discrimination (PSD). Set to `NaN` when PSD simulation is not available (i.e., when the drift-time map or current-pulse model is missing). This is distinct from the usability flags that track LEGEND-200 data quality. |
| `is_single_site`  | `Array` | —     | Boolean PSD flag. `True` if `aoe` falls within the single-site acceptance window defined by cut values extracted from LEGEND-200 data (`psdcuts.aoe.low_side` to `psdcuts.aoe.high_side`).                                                                                                             |

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

| Field            | Type              | Units | Description                                                                                       |
| ---------------- | ----------------- | ----- | ------------------------------------------------------------------------------------------------- |
| `is_good`        | `VectorOfVectors` | —     | Boolean. `True` if the PSD usability flag is valid in LEGEND-200 data. Variable-length per event. |
| `aoe`            | `VectorOfVectors` | —     | A/E classifier values forwarded from the `hit` tier. Variable-length per event.                   |
| `has_aoe`        | `VectorOfVectors` | —     | Boolean. `True` if the A/E value is not `NaN` (i.e. PSD was computed). Variable-length per event. |
| `is_single_site` | `VectorOfVectors` | —     | Boolean PSD flag forwarded from the `hit` tier. Variable-length per event.                        |

The table `geds/psd` contains the single template based PSD simulation, while
`psd_psl` contains the pulse-shape-library based simulation.

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

| Field           | Type     | Units | Description                                                                                                                        |
| --------------- | -------- | ----- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `nr_sim_events` | `Scalar` | —     | Total number of simulated primary events (`primaries_per_job` × `number_of_jobs`). Used to normalise PDFs to physical event rates. |

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
