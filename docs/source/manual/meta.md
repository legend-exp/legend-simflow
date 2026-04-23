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
```

- `add_random_coincidences` (bool) — when `true`, random-coincidence (RC) SiPM
  data (taken from `l200data` evt/pet tiers) is mixed in during event building.
- `geds_energy_thr_kev` (int) — HPGe hit energy threshold in keV; hits below
  this value are discarded.
- `spms_energy_thr_pe` (int) — SiPM hit threshold in photoelectrons; hits below
  this value are discarded.
- `buffer_len` (str) — LH5 read chunk size (e.g. `"50*MB"`). Controls memory
  usage during processing; does not affect the output.

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
```

- `buffer_len` (str) — LH5 read chunk size (e.g. `"500*MB"`). Controls memory
  usage during processing; does not affect the output.

## `pars/` — simulation parameters

Metadata is organized in this directory by experimental configuration (first
level) and detector type (second level), mirroring the `tier/` structure.

### Drift time map settings

A single shared YAML file (applies to all detectors and voltages) that overrides
simulation control parameters for the
[`build_hpge_drift_time_map`](../api/snakemake_rules.md) rule. When absent, the
script uses built-in production defaults.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/dtmap/settings.yaml

grid_size_in_mm: 0.5
ssd_refinement_limits: [0.2, 0.1, 0.05, 0.02]
padding: 3
```

| Key                     | Type          | Default                  | Description                                                                                                                                                                                                                         |
| ----------------------- | ------------- | ------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `grid_size_in_mm`       | float         | `0.5`                    | Drift time map grid spacing in mm. Execution time scales quadratically with `1/grid_size_in_mm`.                                                                                                                                    |
| `ssd_refinement_limits` | list of float | `[0.2, 0.1, 0.05, 0.02]` | SSD adaptive-mesh refinement thresholds. Each entry drives one refinement pass; smaller values give a more accurate electric field at higher cost. **Overly coarse values can prevent full detector depletion — change with care.** |
| `padding`               | int           | `3`                      | Number of pixel layers padded around the drift time map boundary to avoid grid edge effects.                                                                                                                                        |

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
resolution parameters. Follows the same structure and four-case logic as
{ref}`eresmod-metadata-dir`.

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

(psdcuts-metadata-dir)=

### PSD cut defaults

An optional validity-based metadata directory providing HPGe-specific PSD cut
values. Follows the same structure and four-case logic as
{ref}`eresmod-metadata-dir`.

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

An optional validity-based metadata directory providing HPGe-specific current
pulse model parameters. Follows the same structure and four-case logic as
{ref}`eresmod-metadata-dir`, but applied per-detector rather than per-run (one
output file per `(runid, hpge_detector)` pair).

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

- `default` _(optional)_ — current pulse model applied to all HPGe detectors not
  listed explicitly.
- `<detector>` _(optional)_ — per-detector override.

Each entry must contain:

- `current_pulse_pars` — mapping of parameter names to their values for the
  current pulse model (`amax`, `mu`, `sigma`, `tail_fraction`, `tau`,
  `high_tail_fraction`, `high_tau`; the last two default to `0` if omitted)
- `mean_aoe` — mean A/E value
- `current_reso` — current resolution (σ) from the noise-fit

See {ref}`hpge-currmod-extraction` for a description of how these files are used
at runtime.
