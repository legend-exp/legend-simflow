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
production parameter database at processing time.

(hpge-eresmod-extraction)=

### HPGe energy resolution (`extract_hpge_observables_models`)

The `extract_hpge_observables_models` rule produces a per-detector YAML file
mapping each HPGe detector to its energy resolution model. It is a _collection_
step: it gathers what it can from the available sources and writes the result to
disk. Completeness is validated downstream in `build_tier_hit` (see
{ref}`build-tier-hit-hpge`).

The exact sources used depend on what is configured:

1. **`l200data` only, no HPGe-specific overrides** — everything available in
   `l200data` is collected, independent of detector status.

2. **`l200data` + HPGe-specific overrides in {ref}`eresmod-metadata-dir`, no
   `default` key** — everything available in `l200data` is collected as in case
   1; the explicitly listed detectors are then overridden.

3. **`l200data` + HPGe-specific overrides with a `default` key** — the metadata
   takes over entirely: every HPGe detector in the channel map is expanded from
   the `default` (with optional per-detector overrides). `l200data` is not
   consulted. All detectors, including non-ON ones, are present in the output.

4. **HPGe-specific overrides with a `default` key, no `l200data`** — same as
   case 3. `l200data` is not required.

(hpge-aoeresmod-extraction)=

### HPGe A/E resolution (`extract_hpge_observables_models`)

The same rule produces a per-detector YAML file for the A/E resolution model,
following the same four-case logic as {ref}`hpge-eresmod-extraction`. The
metadata directory is described in {ref}`aoeresmod-metadata-dir`.

(hpge-psdcuts-extraction)=

### HPGe PSD cuts (`extract_hpge_observables_models`)

The same rule produces a per-detector YAML file for the PSD cut values,
following the same four-case logic as {ref}`hpge-eresmod-extraction`. The
metadata directory is described in {ref}`psdcuts-metadata-dir`.

(hpge-currmod-extraction)=

### HPGe current pulse model (`extract_current_pulse_model`)

The `extract_current_pulse_model` rule produces one YAML file per
`(runid, hpge_detector)` pair holding the current pulse model parameters for
that detector in that run. Unlike the eresmod/aoeresmod/psdcuts extractions
(which are per-run), this step has per-detector granularity. It follows the same
four-case logic as {ref}`hpge-eresmod-extraction`, but cases 3 and 4 bypass
waveform fitting entirely and do not need `l200data`. The metadata directory is
described in {ref}`currmod-metadata-dir`.

## `opt` — optical hit building

:::{todo}

Document the optical hit building step.

:::

## `hit` — hit tier building

(build-tier-hit-hpge)=

### HPGe observable validation (`build_tier_hit`)

`build_tier_hit` reads the per-detector YAML files produced in the `par` step
and validates that every simulated detector has the parameters it needs. The
hard-error vs. fallback policy differs slightly per observable:

| Observable        | Hard error                                                 | Fallback (+ warning) | Fallback key        |
| ----------------- | ---------------------------------------------------------- | -------------------- | ------------------- |
| Energy resolution | ON detector missing entry                                  | `off`/`ac` detector  | `eresmod_default`   |
| A/E resolution    | ON detector with `psd_usability ≠ "missing"` missing entry | all other cases      | `aoeresmod_default` |
| PSD cuts          | ON detector with `psd_usability ≠ "missing"` missing entry | all other cases      | `psdcuts_default`   |

The fallback keys are read from {ref}`hit-tier-settings`. They are never
triggered when a `default` key is present in the corresponding metadata (cases 3
and 4 of the extraction steps above), because all detectors are already covered
in that case.

:::{note}

An ON detector with `psd_usability = "missing"` explicitly signals that PSD data
are unavailable for that detector (e.g. a known hardware issue). It is therefore
acceptable to fall back to the default A/E resolution and PSD cuts for such
detectors rather than raising a hard error.

:::

## `evt` — event building

:::{todo}

Document the event building step.

:::

## `cvt` — event concatenation

:::{todo}

Document the event concatenation step.

:::

## `pdf` — PDF generation

:::{todo}

Document the PDF generation step.

:::
