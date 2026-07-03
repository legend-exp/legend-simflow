# Data names and paths

## File naming conventions

The file naming convention for simulation output files is:

```
{experiment}-{simid}-job_{jobid}-tier_{tier}.{extension}
```

where each label (in curly brackets `{}`) should be alphanumeric (including
underscores: `_`) when possible. Avoid dashes (`-`) and other special characters
unless explicitly supported for that specific field.

- `experiment` — name or label for the experimental configuration being
  simulated.

- `simid` — stands for "simulation identifier", i.e. a string to uniquely label
  a simulation. Typically specifies the physical process being simulated and the
  experiment's components involved.

- `jobid` — stands for "(simulation) job identifier". It's a zero-padded integer
  that labels independent jobs across which the simulation is split.

- `tier` — the three-character label of the tier. At the moment the simflow
  supports `vtx`, `stp`, `opt`, `hit`, `evt`, `cvt` and `pdf` tiers.

- `extension` — file extension. `lh5` for LEGEND HDF5 files, `gdml` for GDML
  geometry files, `yaml` for plain-text YAML configuration files, `log` for log
  files.

# Metadata

The workflow configuration metadata is stored in
[legend-simflow-config](https://github.com/legend-exp/legend-simflow-config).
Documentation about which metadata is stored there (e.g. which LEGEND
experimental configuration are supported) can be found in the `README.md`.

In this section, the specification of the metadata format is documented.

## `tier/` static tier configuration

Metadata is organized in this directory by tier (first level) and experimental
configuration (second level).

(vtx-tier-meta)=

### `vtx` tier

This section specifies how to configure simulations for the `vtx` tier,
consisting of simulated event vertices for the `stp` tier.

The configuration folder must contain the following files:

- `tier/vtx/{experiment}/simconfig.yaml`

(vtx-simconfig.yaml)=

#### `simconfig.yaml`

This file defines a dictionary of commands used to generate vertices. The keys
in this dictionary can be referenced in the `stp` tier configuration metadata
({ref}`simconfig.yaml`).

The configuration block must contain the `command` key, defining the command
block that should be used to generate the vertices. Then command string must
contain the following variables, that will be automatically substituted by the
simflow at runtime:

- `GDML_FILE`: path to the GDML file defining the geometry that is being
  simulated. This is typically a required input of vertex generators.
- `OUTPUT_FILE`: output file where the vertices will be saved. This file will be
  an input of consumer `stp`-tier jobs.
- `N_EVENTS`: number of vertices to generate.

Example:

```yaml
hpge_surface:
  command: >-
    revertex hpge-surf-pos --detectors [VB]* --surface-type nplus --gdml
    {INPUT_FILE} --out-file {OUTPUT_FILE} --n-events {N_EVENTS}
```

### `stp` tier

This section specifies how to configure
[_remage_](https://remage.readthedocs.io/) simulations for the `stp` tier. The
same conventions also apply to other tiers when relevant.

The configuration folder must contain the following files:

- `tier/stp/{experiment}/simconfig.yaml`
- `tier/stp/{experiment}/generators.yaml`
- `tier/stp/{experiment}/confinement.yaml`

(simconfig.yaml)=

#### `simconfig.yaml`

This is the main configuration file, which defines the set of _remage_
simulations to run. It is a mapping from `simid` to a configuration block, which
configures how to generate the _remage_ macro for that simulation.

:::{important}

`simid` keys must only contain word characters and hyphens, matching the pattern
`[-\w]+` (letters `a–z`, `A–Z`, digits `0–9`, underscores `_`, and hyphens `-`).
In particular, **dots (`.`) are forbidden**: they are the separator between tier
and simid in the `simlist` format (`<tier>.<simid>`), so a dot inside a `simid`
would break parsing.

In this context hyphens are technically allowed by validation, but naming with
**snake case** (letters, digits, underscores only) is still recommended for
clarity, because hyphens are field separators in output file names (e.g.
`{experiment}-{simid}-job_{jobid}-tier_{tier}.lh5`).

:::

Supported fields per `simid`:

- `template` — path to a macro template file. The template must include
  placeholders used by the workflow (see {ref}`macro-templates-subst`). The
  special variable `$_` is substituted with the path to the directory that
  contains the configuration file.
- `generator` — a string:
  - formatted as `~defines:NAME`, where `NAME` is defined in
    {ref}`generators.yaml`.
  - formatted as `~vertices:NAME`, where `NAME` references a vertices simulation
    from the `vtx` tier (see {ref}`vtx-tier-meta`). When vertices are used as
    the generator they carry vertex position _and_ kinematics, so the
    `confinement` key (see below) is forbidden.

- `confinement` — one of:
  - `~defines:NAME` to reference a confinement block in {ref}`confinement.yaml`
  - `~volumes.bulk:PATTERN` to confine to physical volumes matching `PATTERN`
  - `~volumes.surface:PATTERN` to sample on the surface of volumes matching
    `PATTERN`
  - a list of the above strings to combine multiple volume patterns
  - `~vertices:NAME` to used the vertex positions similated by the `vtx` tier
    generator `NAME` (see {ref}`vtx-tier-meta`).
  - `~function:NAME` to use a user-defined function to generate macro commands.
    `NAME` should be in a format:
  ```
    module.function(<...>,*args,**kwargs)
  ```
  see {func}`legendsimflow.commands.get_confinement_from_function` for more
  details. This function should return a list of the _remage_ macro commands.
- `primaries_per_job` — integer, the number of primaries per job; becomes
  `N_EVENTS` in the macro file.
- `number_of_jobs` — integer, how many jobs to split the simulation into.
- `macro_substitutions`: Optional mapping of additional placeholders to values
  to inject into the template (e.g. `HPGE_ENERGY_THRESHOLD: 450 keV`).
- `geom_config_extra`: Optional nested structure to tweak geometry configuration
  for this `simid`. This configuration block is injected unmodified to the
  geometry tooling (currently
  [legend-pygeom-l200](https://legend-pygeom-l200.readthedocs.io)).

Example:

```yaml
hpge_bulk_Rn222_to_Po214:
  template: $_/template.mac
  generator: ~defines:Rn222_to_Po214
  confinement: ~defines:hpge_bulk
  primaries_per_job: 10_000
  number_of_jobs: 4
```

(generators.yaml)=

#### `generators.yaml`

Defines reusable generator command snippets (Geant4 macro lines). Each key can
be:

- A list of strings (recommended), or
- A single string containing multiple lines.

Example block simulating the $^{238}$U decay chain segment from $^{222}$Rn to
$^{214}$Po:

```yaml
Rn222_to_Po214:
  - /RMG/Generator/Select GPS
  - /gps/particle ion
  - /gps/energy 0 eV
  - /gps/ion 86 222
  - /process/had/rdm/nucleusLimits 214 222 82 86
```

In {ref}`simconfig.yaml`, reference this via `generator`:
`~defines:Rn222_to_Po214`.

(confinement.yaml)=

#### `confinement.yaml`

Defines reusable
[_remage_ vertex confinement](https://remage.readthedocs.io/en/stable/manual/confinement.html)
snippets. Each key is a list of lines that will be inserted into the macro when
referenced with `~defines:`.

Example:

```yaml
hpge_bulk:
  - /RMG/Generator/Confine Volume
  - /RMG/Generator/Confinement/Physical/AddVolume V.*
  - /RMG/Generator/Confinement/Physical/AddVolume B.*
```

Alternatively, {ref}`simconfig.yaml` supports direct volume confinement without
an entry in {ref}`confinement.yaml`:

- `~volumes.bulk:REGEX` translates to:

  ```
  /RMG/Generator/Confine Volume
  /RMG/Generator/Confinement/Physical/AddVolume REGEX
  ```

- `~volumes.surface:REGEX` additionally sets:
  ```
  /RMG/Generator/Confinement/SampleOnSurface true
  ```

A list of such tokens combines multiple volume patterns.

(macro-templates-subst)=

#### Macro templates and substitutions

Template macros are standard _remage_/Geant4 macro files that can contain
variable placeholders that the workflow substitutes:

- `$GENERATOR`: Replaced with the generator block content.
- `$CONFINEMENT`: Replaced with the confinement block content.
- `{SEED}`: A 32-bit random integer generated per job.
- `{N_EVENTS}`: The number of events to simulate for the job, taken from
  `primaries_per_job`.
- Additional placeholders may be provided via the {ref}`simconfig.yaml`
  `macro_substitutions` mapping.

Example template snippets:

```
/RMG/Manager/Randomization/Seed {SEED}
...
$GENERATOR
$CONFINEMENT
...
/run/beamOn {N_EVENTS}
```

The workflow renders the template to a concrete macro and writes it to the
canonical input path for the job. It then builds the
[_remage_ CLI](https://remage.readthedocs.io/en/stable/manual/running.html#command-line-options)
either by passing this macro file with `--macro-substitutions` (`SEED` and
`N_EVENTS`), or by inlining the commands directly when using
[_remage_'s "inline" mode](https://remage.readthedocs.io/en/stable/manual/running.html#executing-commands-in-batch-mode).

### Run partitioning

An important post-processing step of the workflow is to fold detector models
with parameters that vary during the livetime of the experiment (across "data
taking runs"), satisfying the following requirements:

- the contribution of each data taking run over the total in terms of livetime
  must be represented in the final simulated event sample;
- the statistical properties of the simulated event sample must be kept intact.

Ultimately, the simulated event sample must be directly comparable to the
observed event sample. The Simflow will partition the total simulated event
sample (across all simulation jobs) according to the livetime fraction of each
run (taken from
[legend-datasets](https://github.com/legend-exp/legend-datasets)) and apply run
parameters to each partition. The partitions will still all live in the same
table in the output file. A new column named `runid` holding the run index of
each event is appended for convenience.

Run partitioning is applied to both the `opt` and `hit` tiers: the LAr response
is processed with the optical parameters and detector usability of the
corresponding run partition, and the HPGe response with the energy resolution
and PSD parameters of the same partition.

The user selects a list of data taking runs ("run list" or "data set") that they
want to simulate. Runs are specified by run-ids ("run identifiers") in the
format:

```
l200-<period>-<run>-<datatype>
```

where

- `period` is `p03`, `p04`, ...
- `run` is `r000`, `r001`, ...
- `datatype` is `phy`, `cal`, `ssc`, ...

In addition, the Simflow gives the possibility to get the runlist from the
`runlists.yaml` database file stored in the
[legend-datasets](https://github.com/legend-exp/legend-datasets) repository by
prefixing the runid string with `~runlists:`, followed by a dot-separated path
to the database entry. For example:

```
~runlists:valid.phy.p04 -> [
  'l200-p04-r000-phy',
  'l200-p04-r001-phy',
  'l200-p04-r002-phy',
  'l200-p04-r003-phy'
]
```

Run lists passed to the Simflow can include both runids and
runlist-file-queries.

The Simflow supports specifying a global runlist in the main configuration file,
under the field `runlist`:

```{code-block} yaml
:caption: simflow-config.yaml

runlist:
  - l200-p03-r000-phy
  - l200-p03-r001-phy
  - ~runlists:valid.phy.p04
```

:::{note}

Per-simid runlist overrides are configured exclusively in the `hit`-tier
`simconfig.yaml`, even though the same partition is reused by the `opt` tier.
There is no separate `opt`-tier `simconfig.yaml` for this purpose.

```{code-block} yaml
:caption: /.../simprod/config/tier/hit/simconfig.yaml

hpge_bulk_Rn222_to_Po214:
  runlist:
    - l200-p03-r003-phy
    - l200-p03-r004-phy
```

:::

(opt-tier-settings)=

### `opt` tier

This section specifies how to configure the `opt` tier, which processes the
liquid argon scintillation response from the `stp` tier.

```{code-block} yaml
:caption: simprod/config/tier/opt/{experiment}/settings.yaml

scintillator_volume_name: liquid_argon
optmap_per_sipm: true
optmap_scaling_factor: 0.3
photoelectron_resolution_sigma: 0.3
time_resolution_in_ns: 16
max_pes_per_hit_per_sipm: 5
max_pes_per_hit_combined: 100
buffer_len: "10*MB"
```

- `scintillator_volume_name` (str) — name of the scintillator volume in the GDML
  geometry used to identify liquid argon energy depositions (e.g.
  `liquid_argon`).
- `optmap_per_sipm` (bool) — when `true`, photoelectrons are sampled per SiPM
  channel using the per-SiPM optical map; when `false`, the combined map across
  all SiPMs is used.
- `optmap_scaling_factor` (float) — factor multiplied to every map value,
  globally scaling the photoelectron detection probability. Maps produced with
  SiPM PDE set to 1 should use a value equal to the true SiPM PDE.
- `photoelectron_resolution_sigma` (float) — single-photoelectron amplitude
  resolution (σ, relative). Applied as Gaussian smearing to each detected
  photoelectron.
- `time_resolution_in_ns` (float) — SiPM time resolution in nanoseconds.
  Photoelectrons within this window are clustered into a single hit.
- `max_pes_per_hit_per_sipm` (int) — maximum number of photoelectrons per hit
  per SiPM channel (used when `optmap_per_sipm: true`). Limits memory and
  processing time.
- `max_pes_per_hit_combined` (int) — maximum number of photoelectrons per hit
  across all SiPMs combined (used when `optmap_per_sipm: false`).
- `buffer_len` (str) — LH5 read chunk size (e.g. `"10*MB"`). Controls memory
  usage during processing; does not affect the output.

(hit-tier-settings)=

### `hit` tier

This section specifies how to configure the post-processing of the
[_remage_](https://remage.readthedocs.io/) simulations from the `stp` tier. When
absent, accessing any field raises an error.

```{code-block} yaml
:caption: simprod/config/tier/hit/{experiment}/settings.yaml

dead_layer_fraction: 0.5
buffer_len: "500*MB"
simulate_psd: True
simulate_psd_with_psl: False

eresmod_default:
  expression: FWHMLinear
  parameters:
    a: 0.5
    b: 0.001

aoeresmod_default:
  expression: SigmaFit
  parameters:
    a: 0.0001
    b: 0
    c: 1

psdcuts_default:
  aoe:
    low_side: -1.5
    high_side: 3
```

- `dead_layer_fraction` (float) — fraction of the dead layer thickness at which
  the linear ramp in charge collection efficiency starts, between `0` (HPGe
  surface, ramp begins immediately) and `1` (ramp begins only at the full charge
  collection depth, i.e. no partial collection).
- `buffer_len` (str) — LH5 read chunk size (e.g. `"500*MB"`). Controls memory
  usage during processing; does not affect the output.
- `eresmod_default` — energy resolution model applied to non-ON detectors. See
  {ref}`build-tier-hit-hpge` for when this fallback is triggered.
- `aoeresmod_default` — A/E resolution model applied to detectors without a
  per-detector entry. See {ref}`build-tier-hit-hpge` for when this fallback is
  triggered.
- `simulate_psd` (bool, default `True`): enable the single-template A/E PSD
  simulation, written to the `geds/psd` subtable of the `hit` tier.
- `simulate_psd_with_psl` (bool, default `False`): enable the
  pulse-shape-library based PSD simulation (see {ref}`hpge-psl-overview`),
  written to the `geds/psd_psl` subtable. The two flags are independent: enable
  either, both, or neither. Setting both to `False` disables the HPGe PSD
  simulation entirely.
- `psdcuts_default` — PSD cut values applied to detectors without a per-detector
  entry. See {ref}`build-tier-hit-hpge` for when this fallback is triggered.

(evt-tier-settings-meta)=

### `evt` tier

This section specifies how to configure the `evt` tier, which builds physics
events from the hit-level data.

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
  data (taken from `l200data` evt/pet tiers) is mixed in during event building.
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

(cvt-tier-settings)=

### `cvt` tier

This section specifies how to configure the `cvt` tier, which concatenates event
files across simulation jobs.

```{code-block} yaml
:caption: simprod/config/tier/cvt/{experiment}/settings.yaml

buffer_len: "500*MB"
```

- `buffer_len` (str) — LH5 read chunk size (e.g. `"500*MB"`). Controls memory
  usage during processing; does not affect the output.

(pdf-tier-settings)=

### `pdf` tier

This section specifies how to configure the `pdf` tier, which produces
probability density functions from the concatenated event files.

```{code-block} yaml
:caption: simprod/config/tier/pdf/{experiment}/settings.yaml

buffer_len: "500*MB"

# optional: split 1-D PDFs by detector group
detector_groups:
  icpc: "V.*"
  bege: "B.*"
```

- `buffer_len` (str) — LH5 read chunk size (e.g. `"500*MB"`). Controls memory
  usage during processing; does not affect the output.
- `detector_groups` (mapping, optional) — maps group names to Python regex
  strings. Each regex is matched against LEGEND-200 detector names using
  `re.fullmatch`, so `"V.*"` selects all detectors whose names start with `V`.
  When this key is absent, only the implicit `all` group is emitted (equivalent
  to `detector_groups: {all: ".*"}`). Specifying `detector_groups` extends the
  output: every named group is produced in addition to `all`, which is always
  emitted regardless of the config. See {ref}`pdf-tier` for the resulting output
  schema.

## `pars/` — simulation parameters

Metadata is organized in this directory by experimental configuration (first
level) and detector type (second level), mirroring the `tier/` structure.

(ssd-settings-meta)=

### Pulse shape simulation settings

A single shared YAML file (applies to all detectors and voltages) that overrides
[`SolidStateDetectors.jl`](https://juliaphysics.github.io/SolidStateDetectors.jl/stable/)
simulation control parameters for the
[`build_hpge_drift_time_map`](../api/snakemake_rules.md) and
[`build_hpge_pulse_shape_library`](../api/snakemake_rules.md) rules. When
absent, the scripts use built-in production defaults.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/ssd/settings.yaml

grid_size_in_mm: 0.5
ssd_refinement_limits: [0.2, 0.1, 0.05, 0.02]
padding: 3
```

| Key                     | Type          | Default                  | Description                                                                                                                                                                                                                                  |
| ----------------------- | ------------- | ------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `grid_size_in_mm`       | float         | `0.5`                    | Simulation grid spacing in mm. Execution time scales quadratically with `1/grid_size_in_mm`. The built-in default differs per rule (0.5 mm for the drift-time map, 5 mm for the pulse-shape library); an explicit value here overrides both. |
| `ssd_refinement_limits` | list of float | `[0.2, 0.1, 0.05, 0.02]` | SSD adaptive-mesh refinement thresholds. Each entry drives one refinement pass; smaller values give a more accurate electric field at higher cost. **Overly coarse values can prevent full detector depletion — change with care.**          |
| `padding`               | int           | `3`                      | Number of pixel layers padded around the simulated map (drift-time map and pulse-shape library) boundary to avoid grid edge effects.                                                                                                         |

:::{tip}

In test or CI environments, setting `grid_size_in_mm: 10.0` reduces the number
of grid points by a factor of ~400 compared to the 0.5 mm production default,
cutting script runtime from many minutes to seconds.

:::

(eresmod-metadata-dir)=

### Energy resolution model defaults

An optional validity-based metadata directory providing HPGe-specific energy
resolution parameters. When present, it can supplement or fully replace
`l200data` as the source of energy resolution parameters — enabling simulations
for experiments that have not yet collected data (e.g. LEGEND-1000). The
structure follows the same validity-based format as
`pars/{experiment}/geds/opv/`.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/eresmod/l200-p03-r%-T%-all-eresmod.yaml

default:
  expression: FWHMLinear
  parameters:
    a: 0.5
    b: 0.001

# optional per-detector override
V02160A:
  expression: FWHMLinear
  parameters:
    a: 0.3
    b: 0.0009
```

- `default` _(optional)_ — energy resolution model applied to all HPGe detectors
  not listed explicitly.
- `<detector>` _(optional)_ — per-detector override; key is the detector name as
  it appears in the channel map (e.g. `V02160A`).

Each entry must contain:

- `expression` — name of the energy resolution function (e.g. `FWHMLinear`)
- `parameters` — mapping of parameter names to their values

See {ref}`hpge-eresmod-extraction` for a description of how these files are used
at runtime.

(aoeresmod-metadata-dir)=

### A/E resolution model defaults

An optional validity-based metadata directory providing HPGe-specific A/E
resolution parameters. Follows the same structure as {ref}`eresmod-metadata-dir`
and the same {ref}`source-resolution logic <par-collection-model>`.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/aoeresmod/l200-p03-r%-T%-all-aoeresmod.yaml

default:
  expression: SigmaFit
  parameters:
    a: 0.0001
    b: 0
    c: 1

# optional per-detector override
V02160A:
  expression: SigmaFit
  parameters:
    a: 0.0002
    b: 0
    c: 1
```

- `default` _(optional)_ — A/E resolution model applied to all HPGe detectors
  not listed explicitly.
- `<detector>` _(optional)_ — per-detector override.

Each entry must contain:

- `expression` — name of the A/E resolution function (e.g. `SigmaFit`)
- `parameters` — mapping of parameter names to their values

See {ref}`hpge-aoeresmod-extraction` for a description of how these files are
used at runtime.

(aoemeanmod-metadata-dir)=

### A/E mean energy-dependence model

A **mandatory** validity-based metadata directory providing the HPGe A/E mean as
a function of energy. Unlike {ref}`eresmod-metadata-dir`,
{ref}`aoeresmod-metadata-dir` and {ref}`psdcuts-metadata-dir`, this model cannot
be extracted from `l200data`, so the metadata is the only source and must be
present for every simulated run. It is read directly when building the `hit`
tier.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/aoemeanmod/l200-p03-r%-T%-all-aoemeanmod.yaml

default:
  single_template:
    expression: a * x + b
    pars:
      a: -0.000001
      b: 1.0
  psl:
    expression: a * x + b
    pars:
      a: -0.000002
      b: 1.0

# optional per-detector override
V02160A:
  single_template:
    expression: a * x + b
    pars:
      a: -0.000003
      b: 1.0
  psl:
    expression: a * x + b
    pars:
      a: -0.000004
      b: 1.0
```

- `default` _(optional)_: model applied to all HPGe detectors not listed
  explicitly. When present it is expanded across all `geds` detectors in the
  channel map, mirroring {ref}`aoeresmod-metadata-dir`.
- `<detector>` _(optional)_: per-detector override.

Each entry must provide a `single_template` and a `psl` block (one per PSD
simulation type), and each block must contain:

- `expression`: a Python expression for the A/E mean as a function of the energy
  `x` (in keV), using the parameters `a` and `b`
- `pars`: mapping of the parameters `a` and `b` to their values

A detector with neither an explicit entry nor a `default` falls back to a flat
A/E mean of 1 (a warning is logged for `on` detectors).

(psdcuts-metadata-dir)=

### PSD cut defaults

An optional validity-based metadata directory providing HPGe-specific PSD cut
values. Follows the same structure as {ref}`eresmod-metadata-dir` and the same
{ref}`source-resolution logic <par-collection-model>`.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/psdcuts/l200-p03-r%-T%-all-psdcuts.yaml

default:
  aoe:
    low_side: -1.5
    high_side: 3.0

# optional per-detector override
V02160A:
  aoe:
    low_side: -1.4
    high_side: 2.9
```

- `default` _(optional)_ — PSD cut values applied to all HPGe detectors not
  listed explicitly.
- `<detector>` _(optional)_ — per-detector override.

Each entry must contain:

- `aoe.low_side` — lower A/E classifier cut (in units of A/E σ)
- `aoe.high_side` — upper A/E classifier cut (in units of A/E σ)

See {ref}`hpge-psdcuts-extraction` for a description of how these files are used
at runtime.

(currmod-metadata-dir)=

### Current pulse model defaults

:::{note}

The current-pulse model drives the single-template A/E PSD simulation, enabled
with `simulate_psd: True` in {ref}`hit-tier-settings` (default `True`). Set it
to `False` to disable the single-template PSD path. See
{ref}`build-tier-hit-hpge`.

:::

An optional validity-based metadata directory providing HPGe-specific current
pulse model parameters. Follows the same structure as
{ref}`eresmod-metadata-dir` and the same
{ref}`source-resolution logic <par-collection-model>`, but applied per-detector
rather than per-run (one output file per `(runid, hpge_detector)` pair).

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/currmod/l200-p03-r%-T%-all-currmod.yaml

default:
  current_pulse_pars:
    amax: 1.0
    mu: 0.0
    sigma: 0.1
    tail_fraction: 0.5
    tau: 0.02
    high_tail_fraction: 0.0
    high_tau: 0.0
  mean_aoe: 1.0
  current_reso: 0.01

# optional per-detector override
V02160A:
  current_pulse_pars:
    amax: 1.0
    mu: 0.0
    sigma: 0.12
    tail_fraction: 0.55
    tau: 0.025
    high_tail_fraction: 0.0
    high_tau: 0.0
  mean_aoe: 0.98
  current_reso: 0.012
```

- `default` _(optional)_ — current-pulse model applied to all HPGe detectors not
  listed explicitly.
- `<detector>` _(optional)_ — per-detector override.

Each entry must contain:

- `current_pulse_pars` — mapping of parameter names to their values for the
  current-pulse model (`amax`, `mu`, `sigma`, `tail_fraction`, `tau`,
  `high_tail_fraction`, `high_tau`; the last two default to `0` if omitted)
- `mean_aoe` — mean A/E value
- `current_reso` — current resolution (σ) from the noise-fit

See {ref}`hpge-currmod-extraction` for a description of how these files are used
at runtime.

(superpulses-settings-meta)=

### Superpulse settings

:::{note}

These settings apply only when the pulse-shape-library (PSL) PSD simulation is
enabled with `simulate_psd_with_psl: True` in {ref}`hit-tier-settings` (default
`False`). See {ref}`hpge-psl-overview`.

:::

A static YAML settings file that controls the superpulse building step
(`build_superpulses_from_data`). When absent, the script uses built-in
production defaults. Settings are read via
`get_par_settings(config, "superpulses")` from
`simprod/config/pars/{experiment}/geds/superpulses/`.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/superpulses/settings.yaml

build_per_runid: false
min_number_wfs: 10
target_wfs: 100
chi2_threshold: 3
t0_field: spms/event_t0
end_time_field: geds/psd/low_aoe/time
drift_time_slices: "1000:200:2000"
evt_tier_name: pet
max_files: null
charge_output: wf_pz_win
curr_output: curr_av
energy_output: cuspEmax
```

| Key                 | Type        | Default                 | Description                                                                                                                                                                                  |
| ------------------- | ----------- | ----------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `build_per_runid`   | bool        | `false`                 | When `true`, produce one superpulse LH5 file per `(runid, detector)` instead of one file per detector accumulating all runs. See {ref}`hpge-superpulses-extraction`.                         |
| `min_number_wfs`    | int         | `10`                    | Minimum waveforms required to form a superpulse for a slice; slices with fewer events are skipped.                                                                                           |
| `target_wfs`        | int         | `100`                   | Target number of waveforms to collect per slice before stopping accumulation.                                                                                                                |
| `chi2_threshold`    | float       | `3`                     | Reduced chi2 threshold for the self-similarity cut; waveforms above this value are rejected before computing the final superpulse.                                                           |
| `t0_field`          | str         | `spms/event_t0`         | Event-level field used as the start time for drift-time calculation. Events where this field is NaN are discarded; for `spms/event_t0` this removes events without a coincident SiPM signal. |
| `end_time_field`    | str         | `geds/psd/low_aoe/time` | Event-level field used as the end time for drift-time calculation.                                                                                                                           |
| `drift_time_slices` | str         | `"1000:200:2000"`       | Drift-time bins as `start:step:stop` in ns; the default creates 200 ns-wide bins from 1000 to 2000 ns.                                                                                       |
| `evt_tier_name`     | str         | `pet`                   | Name of the evt tier to read from the data production (e.g. `pet` or `evt`).                                                                                                                 |
| `max_files`         | int or null | `null`                  | If set, limits the number of raw and evt files processed per run (useful for testing).                                                                                                       |
| `charge_output`     | str         | `wf_pz_win`             | DSP processing chain output name for the charge waveform.                                                                                                                                    |
| `curr_output`       | str         | `curr_av`               | DSP processing chain output name for the current waveform.                                                                                                                                   |
| `energy_output`     | str         | `cuspEmax`              | DSP processing chain output name for the energy estimator used in waveform normalisation.                                                                                                    |

(elecmod-metadata-dir)=

### Electronics-response model defaults

:::{note}

These settings apply only when the pulse-shape-library (PSL) PSD simulation is
enabled with `simulate_psd_with_psl: True` in {ref}`hit-tier-settings` (default
`False`). See {ref}`hpge-psl-overview`.

:::

An optional validity-based metadata directory providing HPGe-specific
electronics-response model parameters. When a `default` key is present, the
data-driven fit in `extract_electronics_model_pars` is bypassed entirely and the
metadata values are written directly to the output YAML. This makes PSL-based
simulations possible without access to LEGEND-200 data (e.g. for LEGEND-1000
studies). The structure follows the same validity-based format as
{ref}`eresmod-metadata-dir`.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/elecmod/l200-p03-r%-T%-all-elecmod.yaml

default:
  sigma: 10.0
  tau: 50.0

# optional per-detector override
V02160A:
  sigma: 12.0
  tau: 45.0
```

- `default` _(optional)_ — electronics-response model parameters applied to all
  HPGe detectors not listed explicitly. When present, also removes the
  requirement for `build_superpulses_from_data` as an input to
  `extract_electronics_model_pars`.
- `<detector>` _(optional)_ — per-detector override; key is the detector name as
  it appears in the channel map (e.g. `V02160A`).

Each entry must contain:

- `sigma` — Gaussian sigma of the digitizer bandwidth in ns.
- `tau` — exponential decay constant of the preamplifier response in ns.

See {ref}`hpge-elecmod-extraction` for a description of how these files are used
at runtime.

(modeling-settings-meta)=

### HPGe modeling settings

An optional single shared YAML file (applies to all detectors) that tunes the
automatic gate deciding which HPGe detectors are eligible for drift-time map and
current-pulse model generation (see {ref}`hpge-modeling-criteria`). When absent,
the built-in defaults are used.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/modeling/settings.yaml

min_voltage_above_depletion_in_V: 100
```

| Key                                | Type | Default | Description                                                                                                                                                    |
| ---------------------------------- | ---- | ------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `min_voltage_above_depletion_in_V` | int  | `100`   | Minimum margin (in V) by which a detector's operational voltage must exceed its depletion voltage to be considered modelable. Detectors below it are excluded. |

(skip-metadata-dir)=

### Manual HPGe skip-list

An optional validity-based metadata directory listing HPGe detectors that should
be excluded from drift-time map and current-pulse model generation, regardless
of their status in the channel map.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/skip/l200-p03-r%-T%-all-skip.yaml

V02160A: "broken impurity profile in legend-metadata"
B00091B: "asymmetric geometry not yet supported"
```

Each entry is a mapping of detector name (as it appears in the channel map, e.g.
`V02160A`) to a free-form reason string describing why the detector is excluded.
The reason is written to the workflow log as a WARNING when the skip is applied.

Detectors listed here are removed from the "modelable" HPGe list for the
matching runs. They will not get a drift-time map nor a current-pulse model
produced. The validity rules are the same as those of the other `geds/`
parameter directories (see {ref}`eresmod-metadata-dir`).

A detector that is manually skipped is treated identically to one that lacks a
drift-time map or current-pulse model for any other reason: PSD output columns
are filled with NaN and the fallback A/E resolution and PSD cuts
(`aoeresmod_default` / `psdcuts_default`) are used. No hard error is raised. See
{ref}`build-tier-hit-hpge` for the full fallback policy.
