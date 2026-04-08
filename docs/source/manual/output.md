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

| Field             | Type    | Units | Description                                                                                                                                                                                                                                                                         |
| ----------------- | ------- | ----- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `energy`          | `Array` | keV   | Reconstructed energy after smearing with the detector energy resolution. Computed from the sum of active energy depositions (weighted by the dead-layer activeness model).                                                                                                          |
| `drift_time_amax` | `Array` | ns    | Drift time at the maximum-current (A) position of the simulated current pulse. Set to `NaN` when no drift-time map or current-pulse model is available for the detector.                                                                                                            |
| `aoe_raw`         | `Array` | —     | Raw A/E value: maximum current amplitude divided by energy. The maximum current (A) is obtained from the simulated current pulse, constructed from individual hit drift times and a current-pulse model, and includes electronic noise effects.                                     |
| `aoe`             | `Array` | —     | A/E classifier value: `(aoe_raw - 1) / aoe_resolution`. Used for pulse-shape discrimination (PSD). Set to `NaN` when PSD simulation is not available (drift-time map or current-pulse model missing); this is distinct from the usability flags that track LEGEND-200 data quality. |
| `is_single_site`  | `Array` | —     | Boolean PSD flag. `True` if `aoe` falls within the single-site acceptance window defined by cut values extracted from LEGEND-200 data (`psdcuts.aoe.low_side` to `psdcuts.aoe.high_side`).                                                                                          |
| `period`          | `Array` | —     | Data-taking period number extracted from the run identifier (numeric encoding).                                                                                                                                                                                                     |
| `run`             | `Array` | —     | Data-taking run number extracted from the run identifier (numeric encoding).                                                                                                                                                                                                        |
| `usability`       | `Array` | —     | Encoded detector usability status for this run (e.g. `on`, `off`, `ac`). Decode with {func}`legendsimflow.metadata.decode_usability`. See the detector status flags in `legend-metadata/datasets/statuses`.                                                                         |
| `psd_usability`   | `Array` | —     | Encoded PSD usability flag (e.g. `valid`). Indicates whether PSD parameters are valid in LEGEND-200 data for this detector and run. Decode with {func}`legendsimflow.metadata.decode_psd_usability`.                                                                                |

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
| `aoe`             | `VectorOfVectors` | —     | A/E classifier values forwarded from the `hit` tier. Variable-length per event.                                         |
| `has_aoe`         | `VectorOfVectors` | —     | Boolean. `True` if the A/E value is not `NaN` (i.e. PSD was computed). Variable-length per event.                       |
| `is_single_site`  | `VectorOfVectors` | —     | Boolean PSD flag forwarded from the `hit` tier. Variable-length per event.                                              |
| `multiplicity`    | `Array`           | —     | Number of HPGe hits above threshold per event. Scalar per event.                                                        |

#### `geds/psd/` — PSD quality

| Field     | Type              | Units | Description                                                                                       |
| --------- | ----------------- | ----- | ------------------------------------------------------------------------------------------------- |
| `is_good` | `VectorOfVectors` | —     | Boolean. `True` if the PSD usability flag is valid in LEGEND-200 data. Variable-length per event. |

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
