# Pipeline steps

This page describes what each step of the Simflow pipeline does and how it uses
the metadata and configuration described in [](meta.md).

## `vtx` — event vertex generation

:::{todo} Document the vertex generation step. :::

## `stp` — particle simulation

:::{todo} Document the remage simulation step. :::

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
{ref}`build-tier-hit-energy-resolution`).

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

## `opt` — optical hit building

:::{todo} Document the optical hit building step. :::

## `hit` — hit tier building

(build-tier-hit-energy-resolution)=

### Energy resolution (`build_tier_hit`)

`build_tier_hit` applies energy resolution smearing to each simulated detector.
It reads the per-detector file produced in the `par` step and validates that all
needed detectors are covered:

- **ON detectors** must always have a curve — a missing entry raises a hard
  error.
- **`off` and `ac` detectors** fall back to `eresmod_default` from
  {ref}`hit-tier-settings` with a warning. This fallback is never triggered when
  a `default` key is present in the eresmod metadata (cases 3 and 4 above),
  because all detectors are already covered.

## `evt` — event building

:::{todo} Document the event building step. :::

## `cvt` — event concatenation

:::{todo} Document the event concatenation step. :::

## `pdf` — PDF generation

:::{todo} Document the PDF generation step. :::
