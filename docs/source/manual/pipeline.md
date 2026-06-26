# Pipeline steps

This page describes what each step of the Simflow pipeline does and how it uses
the metadata and configuration described in [](meta.md).

## `vtx` — event vertex generation

:::{todo}

Document the vertex generation step.

:::

## `stp` — particle simulation

:::{todo}

Document the remage simulation step.

:::

## `par` — parameter extraction

The `par` step extracts per-detector observable models for each data taking run
before the `hit` tier is built. The results are written to per-run YAML files on
disk so that `build_tier_hit` does not need to load the (potentially large) data
production parameter database at processing time. The on-disk locations of the
files produced by this step are listed in {ref}`par-output-files`.

(par-collection-model)=

### Parameter source resolution

Several `par`-tier rules build their per-detector output by combining two
possible sources: the LEGEND-200 data production (`l200data`) and an optional
HPGe-specific metadata override directory (a validity-based directory, one per
observable). Each such rule is a _collection_ step, and which source is used
depends on what is configured:

1. **`l200data` only, no HPGe-specific overrides** — everything available in
   `l200data` is collected, independent of detector status.

2. **`l200data` + HPGe-specific overrides without a `default` key** — everything
   available in `l200data` is collected as in case 1; the explicitly listed
   detectors are then overridden.

3. **`l200data` + HPGe-specific overrides with a `default` key** — the metadata
   takes over entirely: every HPGe detector in the channel map is expanded from
   the `default` (with optional per-detector overrides). `l200data` is not
   consulted. All detectors, including non-ON ones, are present in the output.

4. **HPGe-specific overrides with a `default` key, no `l200data`** — same as
   case 3. `l200data` is not required.

Completeness is validated downstream in `build_tier_hit` (see
{ref}`build-tier-hit-hpge`). Cases 3 and 4 cover every detector, so the
`build_tier_hit` fallback defaults are never triggered.

(hpge-eresmod-extraction)=

### HPGe energy resolution

The `extract_hpge_observables_models` rule produces a per-detector YAML file
mapping each HPGe detector to its energy resolution model, following the
{ref}`source-resolution logic <par-collection-model>`. Its HPGe-specific
metadata override directory is described in {ref}`eresmod-metadata-dir`.

(hpge-aoeresmod-extraction)=

### HPGe A/E resolution

The same rule produces a per-detector YAML file for the A/E resolution model,
following the {ref}`source-resolution logic <par-collection-model>`. The
metadata directory is described in {ref}`aoeresmod-metadata-dir`.

(hpge-psdcuts-extraction)=

### HPGe PSD cuts

The same rule produces a per-detector YAML file for the PSD cut values,
following the {ref}`source-resolution logic <par-collection-model>`. The
metadata directory is described in {ref}`psdcuts-metadata-dir`.

(hpge-currmod-extraction)=

### HPGe current-pulse model

The `extract_current_pulse_model` rule produces one YAML file per
`(runid, hpge_detector)` pair holding the current-pulse model parameters for
that detector in that run. It follows the
{ref}`source-resolution logic <par-collection-model>`, with two differences: it
has per-detector rather than per-run granularity, and cases 3 and 4 bypass the
waveform fitting entirely (so they do not need `l200data`). The metadata
directory is described in {ref}`currmod-metadata-dir`.

### Manual HPGe skip-list

In addition to the automatic exclusions performed by
`gen_list_of_hpges_valid_for_modeling` (usability `on` + valid crystal metadata
required), individual detectors can be removed from the modelable HPGe list for
a given set of runs via the {ref}`skip-metadata-dir` directory. A WARNING is
emitted for each skipped detector, including the reason string and the run
identifier.

(hpge-psl-overview)=

### Pulse-shape library (PSL) chain

:::{note}

This chain is built only when `simulate_psd_with_psl: True` is set in
{ref}`hit-tier-settings` (default `False`). With the flag disabled the realistic
PSL is not requested by `build_tier_hit` and none of the `par`-tier rules in
this chain run.

:::

The PSL chain produces the realistic per-detector waveform library consumed by
`build_tier_hit` for PSL-based pulse-shape discrimination. It proceeds through
four stages, each a separate `par`-tier rule:

1. Compute the ideal pulse-shape libraries with
   [`SolidStateDetectors.jl`](https://juliaphysics.github.io/SolidStateDetectors.jl/stable/)
   (simulation-derived, no real data needed).
2. Build data superpulses (average charge and current waveforms per drift-time
   slice) from LEGEND-200 data.
3. Fit the electronics transfer function (Gaussian digitizer bandwidth $\sigma$
   and preamplifier decay $\tau$) by comparing processed ideal PSL waveforms to
   the data superpulses.
4. Convolve the ideal PSL with the fitted transfer function to produce the
   realistic PSL, which also carries the per-pixel drift-time map used for
   PSL-based PSD.

The realistic PSL is then consumed by `build_tier_hit` as the
pulse-shape-library input to reboost's HPGe PSD routines
({func}`reboost.hpge.utils.get_hpge_pulse_shape_library`).

(hpge-dtmap-extraction)=

### HPGe drift-time maps

:::{note}

The standalone drift-time map feeds the **single-template** PSD path
(`simulate_psd`); the PSL path derives its per-pixel drift times from the
realistic PSL instead. It is documented here because it is, like the ideal PSL,
an `SolidStateDetectors.jl` simulation sharing the same
{ref}`ssd-settings-meta`.

:::

`build_hpge_drift_time_map` runs a Julia script backed by
`SolidStateDetectors.jl` for each `(hpge_detector, hpge_voltage)` pair. The
script reads crystal geometry from `legend-metadata` (diode and crystal YAML
files), solves the electric field, and records the drift time from each $(r, z)$
grid point to the readout contact. The output is one LH5 file per
`(detector, voltage)` pair, in the $(r, z)$-field format consumed by reboost
({func}`reboost.hpge.psd.drift_time`). See the
[`build_hpge_drift_time_map`](../api/snakemake_rules.md) rule reference for the
output fields.

`merge_hpge_drift_time_maps` copies the top-level LH5 objects from all
per-detector files that apply to a given `runid` into a single merged LH5 file
keyed by detector name. If no drift-time maps exist for a run, an empty HDF5
file is created.

`plot_hpge_drift_time_maps` produces a validation PDF for each
`(detector, voltage)` pair from the single-detector output.

The `SolidStateDetectors.jl` (SSD) simulation parameters (grid size, refinement
thresholds, padding) are controlled by {ref}`ssd-settings-meta`.

(hpge-ideal-psl-extraction)=

### Ideal HPGe pulse-shape library

`build_hpge_pulse_shape_library` runs a Julia script backed by
`SolidStateDetectors.jl` for each `(hpge_detector, hpge_voltage)` pair, using
the same crystal geometry as the drift-time map. It simulates the full
charge-collection transient at every $(r, z, \text{angle})$ grid point. These
are un-convolved charge waveforms; the electronics transfer function is applied
in the realistic PSL step. See the
[`build_hpge_pulse_shape_library`](../api/snakemake_rules.md) rule reference for
the output fields.

The SSD simulation parameters are shared with the drift-time map and controlled
by {ref}`ssd-settings-meta`.

(hpge-superpulses-extraction)=

### Data superpulses

`build_superpulses_from_data` builds average ("super") charge and current
waveforms per drift-time slice from LEGEND-200 data for each HPGe detector.
These represent the measured pulse shape at known drift times and are used to
constrain the electronics-response model fit.

**Input data.** For each configured `runid`, the data source is resolved as
follows:

- **Physics runs (`phy`)**: the reference calibration run is looked up and its
  `raw`, `evt` (or `pet`), and channel-map data are used in place of the
  physics-run data. If several physics runs resolve to the same reference
  calibration run, the input files are de-duplicated so that no event is counted
  twice.
- **Non-physics runs** (calibration `cal`, SSC `ssc`, ...): the run's own `raw`
  and `evt` tier files are used directly.

In both cases the script reads `raw`-tier LH5 waveforms and applies the
production DSP configuration.

**Analysis.** Events in a fixed energy range are read from the configured event
tier (`evt_tier_name`, e.g. `evt` or `pet`), passed through standard quality
cuts, and binned by drift time. The drift time is
$t_\text{drift} = t_\text{end} - t_0$, where $t_\text{end}$ and $t_0$ are the
`end_time_field` and `t0_field` configured in {ref}`superpulses-settings-meta`
and events with an undefined $t_0$ are discarded. Within each drift-time slice
(field `drift_time_slices`) the surviving waveforms are averaged into a
preliminary superpulse, a self-similarity $\chi^2$ cut (`chi2_threshold`)
removes outliers, and the final superpulse is the mean of the remaining
waveforms. Slices with fewer than `min_number_wfs` surviving waveforms are
skipped.

See the [`build_superpulses_from_data`](../api/snakemake_rules.md) rule
reference for the per-slice fields. The `build_per_runid` flag (see
{ref}`superpulses-settings-meta`) selects whether superpulses are accumulated
per detector (default) or split per run; the consumer rule
`extract_electronics_model_pars` adjusts its input accordingly.

(hpge-elecmod-extraction)=

### HPGe electronics-response model

`extract_electronics_model_pars` fits the electronics transfer function for each
`(runid, hpge_detector)` pair. The response kernel is the convolution of a
Gaussian (digitizer bandwidth $\sigma$) and a causal exponential decay
(preamplifier response, decay constant $\tau$),

$$
h = h_\text{digi} * h_\text{preamp}, \qquad
h_\text{digi}(t) \propto \exp\!\left[-\frac{1}{2}\left(\frac{t}{\sigma}\right)^2\right], \qquad
h_\text{preamp}(t) \propto e^{-t/\tau}\,\Theta(t),
$$

with $\sigma$ and $\tau$ in ns, $\Theta$ the Heaviside step, and each factor
normalised to unit sum. Convolving an ideal PSL charge waveform with $h$ and
differentiating yields a simulated current waveform. The fit minimises the mean
data-amplitude-weighted RMS between these simulated waveforms and the measured
data superpulses across drift-time slices, using Minuit (MIGRAD).

For each candidate `(sigma, tau)`, ideal PSL waveforms whose drift time falls in
a matching slice are convolved with the electronics kernel, differentiated, and
moving-window-averaged to produce simulated current waveforms. Only slices whose
drift-time centre falls within `dt_range_tuning` are used (defaults to 600-3000
ns), and at most `max_num_superpulses` slices (default `5`, sorted by decreasing
drift time) are selected.

If a `default` key is present in the `elecmod` validity metadata (see
{ref}`elecmod-metadata-dir`), the fitting step is bypassed entirely and the
metadata values are written directly to the output YAML. This is the analogue of
cases 3 and 4 of the {ref}`source-resolution logic <par-collection-model>`; the
difference is that the non-`default` path here is a data-driven fit rather than
a collection from `l200data`.

The fitted parameters (`sigma`, `tau`, and fit diagnostics) are written to a
temporary per-detector YAML; see the
[`extract_electronics_model_pars`](../api/snakemake_rules.md) rule reference for
the full list of keys.

`merge_electronics_model_pars` collects the per-detector YAML files and writes a
single YAML keyed by detector name for each `runid`. This merged file is the
direct input for the realistic PSL convolution step; `sigma` and `tau` are the
only keys required by the consumer.

(hpge-realistic-psl-extraction)=

### Realistic HPGe pulse-shape library

`convolve_hpge_ideal_pulse_shape_lib` selects the ideal PSL for the detector's
operational voltage in that run, reads the fitted `sigma` and `tau` from the
merged electronics model file, and applies the following processing chain:

1. Convolve the ideal charge waveforms with the electronics response kernel
   (Gaussian bandwidth combined with preamplifier decay).
2. Differentiate to obtain current waveforms, then apply a three-stage
   moving-window average.
3. Align waveforms by shifting the current peak to a fixed sample index.
4. Compute the drift time for each (r, z) pixel from the current peak position.
5. Normalise the current amplitudes by the mode of the A/E distribution.

The output for a single `(runid, detector)` pair is in the pulse-shape-library
format consumed by reboost
({func}`reboost.hpge.utils.get_hpge_pulse_shape_library`). See the
[`convolve_hpge_ideal_pulse_shape_lib`](../api/snakemake_rules.md) rule
reference for the output fields. Validation plots (R and Z waveform scans and an
A/E R/Z heatmap) are produced alongside.

`merge_hpge_realistic_psls` copies the per-detector realistic PSL objects into a
single merged LH5 file keyed by detector name for each `runid`. This merged file
is the final output consumed by `build_tier_hit` when `simulate_psd_with_psl` is
enabled.

## `opt` — optical hit building

:::{todo}

Document the optical hit building step.

:::

## `hit` — hit tier building

(build-tier-hit-hpge)=

### HPGe observable validation

:::{note}

HPGe PSD simulation runs in two independent, optional modes selected in
{ref}`hit-tier-settings`: set `simulate_psd: True` (default) for the
single-template A/E estimate (`geds/psd` output) and
`simulate_psd_with_psl: True` for the pulse-shape-library based estimate
(`geds/psd_psl` output, see {ref}`hpge-psl-overview`). Either, both, or neither
can be enabled; with both `False` no PSD columns are written.

:::

`build_tier_hit` reads the per-detector YAML files produced in the `par` step
and validates that every simulated detector has the parameters it needs.

The hard-error vs. fallback policy below applies whenever PSD is simulated and
differs slightly per observable:

| Observable        | Hard error                                                                                             | Fallback (+ warning) | Fallback key        |
| ----------------- | ------------------------------------------------------------------------------------------------------ | -------------------- | ------------------- |
| Energy resolution | ON detector missing entry                                                                              | `off`/`ac` detector  | `eresmod_default`   |
| A/E resolution    | ON detector with `psd_usability ≠ "missing"` and PSD modelable (dtmap + currmod present) missing entry | all other cases      | `aoeresmod_default` |
| PSD cuts          | ON detector with `psd_usability ≠ "missing"` and PSD modelable (dtmap + currmod present) missing entry | all other cases      | `psdcuts_default`   |

The fallback keys are read from {ref}`hit-tier-settings`. They are never
triggered when a `default` key is present in the corresponding metadata (cases 3
and 4 of the {ref}`source-resolution logic <par-collection-model>`), because all
detectors are already covered in that case.

:::{note}

An ON detector with `psd_usability = "missing"` explicitly signals that PSD data
are unavailable for that detector (e.g. a known hardware issue). Similarly, a
detector missing a drift-time map or current-pulse model cannot have its PSD
response simulated at all — the PSD output columns are filled with NaN in that
case. Both situations fall back to the default A/E resolution and PSD cuts
rather than raising a hard error. A detector listed in the
{ref}`skip-metadata-dir` directory is treated identically: it has no drift-time
map or current-pulse model and therefore follows the same NaN PSD output /
fallback A/E resolution and PSD cuts path, with no hard error.

:::

## `evt` — event building

:::{todo}

Document the event building step.

:::

(evt-tier-settings)=

### `simprod/config/tier/evt/{experiment}/settings.yaml` — evt tier settings

A static YAML file with evt-tier-specific settings that apply to all simulations
for a given experiment configuration.

```{code-block} yaml
:caption: simprod/config/tier/evt/{experiment}/settings.yaml

add_random_coincidences: false
geds_energy_thr_kev: 25
spms_energy_thr_pe: 0
buffer_len: "50*MB"
skip_opt: false
skip_hit: false
```

- `add_random_coincidences` (bool) — when `true`, random-coincidence (RC) SiPM
  data is mixed in during event building.
- `geds_energy_thr_kev` (int) — HPGe hit energy threshold in keV; hits below
  this value are discarded.
- `spms_energy_thr_pe` (int) — SiPM hit threshold in photoelectrons; hits below
  this value are discarded.
- `buffer_len` (str) — LH5 read chunk size (e.g. `"50*MB"`). Controls memory
  usage during processing; does not affect the output.
- `skip_opt` (bool, default `false`) — when `true`, the `opt` (SiPM/LAr) tier is
  skipped: the opt Snakemake rule is not run and the evt output contains only
  HPGe data (no `spms` or `coincident/spms` tables).
- `skip_hit` (bool, default `false`) — when `true`, the `hit` (HPGe) tier is
  skipped: the hit Snakemake rule is not run and the evt output contains only
  SiPM data (no `geds` or `coincident/geds` tables).

:::{note}

Setting both `skip_opt` and `skip_hit` to `true` simultaneously is an error.

:::

## `cvt` — event concatenation

:::{todo}

Document the event concatenation step.

:::

## `pdf` — PDF generation

:::{todo}

Document the PDF generation step.

:::
