# Data names and paths

(file-naming)=

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

In this section, the specification of the metadata format is documented. The
defaults of the optional fields of the settings files are listed in
{mod}`legendsimflow.defaults`.

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

Every block must also set the `product` key, which states what the generator
writes. The `stp` tier uses it to choose the _remage_ macro commands for a
`~vertices:` reference (see the `generator` and `confinement` fields in
{ref}`simconfig.yaml`). Allowed values:

- `positions`: vertex positions only, in a `vtx/pos` table. Usable only as an
  `stp`-tier `confinement`.
- `kinematics`: primary kinematics only, in a `vtx/kin` table without positions.
  Usable as an `stp`-tier `generator`, together with a `confinement` that
  supplies the positions.
- `kinematics_and_positions`: primary kinematics and positions in the same
  `vtx/kin` table. Usable as an `stp`-tier `generator` without any
  `confinement`.

Example:

```yaml
hpge_surface:
  command: >-
    revertex hpge-surf-pos --detectors [VB]* --surface-type nplus --gdml
    {INPUT_FILE} --out-file {OUTPUT_FILE} --n-events {N_EVENTS}
  product: positions
```

### `stp` tier

This section specifies how to configure
[_remage_](https://remage.readthedocs.io/) simulations for the `stp` tier. The
same conventions also apply to other tiers when relevant.

The configuration folder must contain the following files:

- `tier/stp/{experiment}/simconfig.yaml`
- `tier/stp/{experiment}/generators.yaml`
- `tier/stp/{experiment}/confinement.yaml`
- `geom/{experiment}-geom-config.yaml`

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
    from the `vtx` tier (see {ref}`vtx-tier-meta`). `NAME` must produce
    kinematics, i.e. its `product` must be `kinematics` or
    `kinematics_and_positions`:
    - with `kinematics_and_positions`, positions are read from the same file, so
      the `confinement` key (see below) is forbidden;
    - with `kinematics`, the file carries no positions, so a `confinement` (of
      any kind) **must** be configured to supply them.

- `confinement` — one of:
  - `~defines:NAME` to reference a confinement block in {ref}`confinement.yaml`
  - `~volumes.bulk:PATTERN` to confine to physical volumes matching `PATTERN`
  - `~volumes.surface:PATTERN` to sample on the surface of volumes matching
    `PATTERN`
  - a list of the above strings to combine multiple volume patterns
  - `~vertices:NAME` to use the vertex positions simulated by the `vtx` tier
    generator `NAME` (see {ref}`vtx-tier-meta`). `NAME` must produce positions,
    i.e. its `product` must be `positions`. This cannot be combined with a
    `~vertices:` generator.
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
  geometry tooling.

Example:

```yaml
hpge_bulk_Rn222_to_Po214:
  template: $_/template.mac
  generator: ~defines:Rn222_to_Po214
  confinement: ~defines:hpge_bulk
  primaries_per_job: 10_000
  number_of_jobs: 4
```

(geom-config.yaml)=

#### `{experiment}-geom-config.yaml`

The geometry configuration file for the experiment. It lives in `geom/`, not in
`tier/stp/`, because the same geometry is used by all tiers. The Simflow copies
it for each `simid`, merges `geom_config_extra` into the copy, and gives the
result to the geometry generator.

```yaml
executable: legend-pygeom-l200
public_geom: true
```

- `executable` (str) — the geometry generator to run, either
  `legend-pygeom-l200` or `legend-pygeom-l1000`. This key is mandatory. The
  generators ignore it. Only simflow reads it.

All other keys are passed to the generator without a change. See the
[legend-pygeom-l200](https://legend-pygeom-l200.readthedocs.io) and
[legend-pygeom-l1000](https://legend-pygeom-l1000.readthedocs.io) manuals for
the keys each one accepts.

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

(run-partitioning)=

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
  geometry used to identify liquid argon energy depositions (e.g. `liquid_argon`
  for LEGEND-200, `liquid_argon_underground` for LEGEND-1000). Not needed if
  `scintillators` is set.
- `scintillators` (dict, optional, default `null`) — the light sources used with
  `light_source: optmap`, `<source name> -> <settings>`. The light of all
  sources reaches the same SiPMs. Each source has:
  - `table` (str) — name of the source's table in the `stp` files. This is the
    name of the scintillator volume registered in the GDML geometry, or the
    ntuple name shared by several registered volumes: e.g. `pen` for all PEN
    plates, which share one optical map. Each source has its own table.
  - `material` (str) — scintillation model, i.e. light yield and time profile,
    as known to reboost: `lar` (liquid argon), `pen` or `fiber` (the core of the
    scintillating fibers, with the time profile of their wavelength shifter).
  - `optmap_scaling_factor` (float or dict, optional) — replaces the top-level
    `optmap_scaling_factor` for this source.
  - `model` (dict, optional) — a light model that replaces the optical map of
    this source, see below. Without it, the source uses its optical map.

  The optical map of source `<name>` is read from `paths.optical_maps.<name>`.
  Without this block, the only source is `lar`, with table
  `scintillator_volume_name`:

  ```{code-block} yaml
  scintillators:
    lar:
      table: liquid_argon
      material: lar
    pen:
      table: pen
      material: pen
      optmap_scaling_factor: 0.5
  ```

  A light model is a Python function that replaces the optical map of a source:
  it gives the probability that a scintillation photon emitted in an energy
  deposition (step) produces a photoelectron. The number of emitted photons
  still comes from the light yield of `material`, as with the optical map. The
  source needs no `paths.optical_maps.<name>`. `model` has two fields:
  - `function` (str) — import path of the function, `<module>.<function>`. The
    module must be importable in the Simflow environment.
  - `kwargs` (dict, optional) — keyword arguments passed to the function, e.g.
    the parameters of the model.

  ```{code-block} yaml
  scintillators:
    lar:
      table: liquid_argon
      material: lar
    pen:
      table: pen
      material: pen
      model:
        function: mypackage.light.pen_detection_probability
        kwargs:
          prob: 0.002
  ```

  The `opt` tier calls the function once per chunk of steps, as
  `f(xloc, yloc, zloc, **kwargs)`. `xloc`, `yloc`, `zloc` hold the position of
  each step in mm, in the global frame of the geometry, as one-dimensional NumPy
  arrays: the function must act on whole arrays. A function with a `geometry`
  keyword argument also gets the simulated GDML geometry, as a
  `pyg4ometry.geant4.Registry`. The function returns an array with one
  probability per step, between 0 and 1, at unit channel efficiency: the
  probability of a photoelectron in any SiPM, like the values of the optical map
  of all SiPMs together. As with the map, the `opt` tier draws the number of
  photoelectrons of each step from a Poisson distribution, with mean the number
  of emitted photons times the probability times `optmap_scaling_factor`. The
  photoelectron times follow the scintillation time profile of `material`. A
  light model gives no light per SiPM channel: with `optmap_per_sipm: true`, the
  `opt` tier leaves out the sources with a model, with a warning. For example:

  ```{code-block} python
  :caption: mypackage/light.py

  import numpy as np


  def pen_detection_probability(xloc, yloc, zloc, *, prob):
      # same probability everywhere
      return np.full_like(xloc, prob)
  ```

  {mod}`legendsimflow.opt_models` provides light models: e.g.
  `legendsimflow.opt_models.uniform` (same probability everywhere) for the PEN
  plates and `legendsimflow.opt_models.l1000_fibers` for the LEGEND-1000 fibers.

  The `opt` tier pools the photoelectrons of all sources that belong to the same
  `stp` event window: hits of the same Geant4 event closer than 10 µs in time,
  as grouped by the `stp` time-coincidence map (TCM). Saturation
  (`max_pes_per_hit_*`), amplitude smearing and time clustering then act on the
  pooled photoelectrons of each SiPM, so light from liquid argon and PEN that
  arrives within `time_resolution_in_ns` forms one pulse. Each window gives one
  output row per SiPM table.

- `light_source` (str, optional, default `optmap`) — where the detected photons
  come from:
  - `optmap`: photoelectrons are sampled from the optical map of each light
    source (`paths.optical_maps.<source>`, see `scintillators`), or from its
    light model, given the energy depositions in its `stp` table.
  - `tracked_photons`: photons are tracked in the simulation, and the photons
    recorded in the `stp` SiPM tables are used. The optical map is not read.
    `max_pes_per_hit_per_sipm` caps the photoelectrons of each SiPM hit, keeping
    the first ones in time. As in the `hit` tier, each row of a SiPM table gives
    one output row, so light that does not come from liquid argon energy
    depositions is kept as well. The output is always per SiPM;
    `optmap_per_sipm` is ignored. Every recorded photon is a photoelectron, so
    the SiPM photon detection efficiency (PDE) must be applied while tracking,
    in the geometry: the PDE curve (`sipm_use_pde_curve`) and the per-channel
    `sipm_efficiencies` of legend-pygeom-l200, set in `geom_config_extra`.
    `optmap_scaling_factor` and `store_expected_pes` are ignored. The `stp` tier
    must be simulated with optical physics and SiPM detectors registered as
    optical detectors. No liquid argon (scintillator) table is needed in the
    `stp` files, and none is read if present.
- `optmap_per_sipm` (bool) — when `true`, photoelectrons are sampled per SiPM
  channel using the per-SiPM optical map (or light model); when `false`, the
  combined map across all SiPMs is used.
- `optmap_scaling_factor` (float or dict, optional if every source in
  `scintillators` sets its own) — factor multiplied to every map value, scaling
  the photoelectron detection probability. The probability given by a light
  model is multiplied by it too. Maps produced with SiPM PDE set to 1 should use
  a value equal to the true SiPM PDE. A scalar is applied to every channel. A
  mapping `<sipm name> -> <factor>` applies a per-channel efficiency instead,
  with a reserved `default` key used for channels absent from it:

  ```{code-block} yaml
  optmap_scaling_factor:
    default: 0.3
    S002: 0.384
    S003: 0.461
  ```

  Channels absent from the mapping fall back to `default`; if neither the
  channel nor `default` is present, processing fails with a `KeyError` naming
  the missing channel.

  Ignored with `light_source: tracked_photons`.

  Per-channel values are equivalent to the `sipm_efficiencies` option of
  [legend-pygeom-l200](https://legend-pygeom-l200.readthedocs.io), which only
  takes effect while photons are tracked and so cannot be used with a pre-built
  optical map.

- `photoelectron_resolution_sigma` (float) — single-photoelectron amplitude
  resolution (σ, relative). Applied as Gaussian smearing to each detected
  photoelectron.
- `time_resolution_in_ns` (float) — SiPM time resolution in nanoseconds.
  Photoelectrons within this window are clustered into a single hit.
- `max_pes_per_hit_per_sipm` (int) — maximum number of photoelectrons per hit
  per SiPM channel (used when `optmap_per_sipm: true` or with
  `light_source: tracked_photons`; `-1` disables it). Limits memory and
  processing time.
- `max_pes_per_hit_combined` (int) — maximum number of photoelectrons per hit
  across all SiPMs combined (used when `optmap_per_sipm: false`).
- `store_expected_pes` (bool, optional, default `false`) — when `true`, also
  store the `expected_pes` field in the `opt` output (see {ref}`opt-tier`) and
  the `spms/expected_pes` field in the `evt` output. Ignored with
  `light_source: tracked_photons`.
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
two_pass_aoe_correction: False
tune_hpge_impurities_on_data: False

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

aoemeanmod_default:
  single_template:
    expression: a * x + b
    pars:
      a: 0
      b: 1
  psl:
    expression: a * x + b
    pars:
      a: 0
      b: 1
```

- `dead_layer_fraction` (float) — fraction of the dead layer thickness at which
  the linear ramp in charge collection efficiency starts, between `0` (HPGe
  surface, ramp begins immediately) and `1` (ramp begins only at the full charge
  collection depth, i.e. no partial collection).
- `store_surface_steps` (bool, default `false`) — store the energy depositions
  near the n+ surface, so that the hit energies can be recomputed for other
  dead-layer models (see {ref}`hit-tier-output`). The `evt` tier still drops the
  hits below its energy threshold, so hits that pass the threshold only for
  another dead-layer model are lost (see {ref}`evt-tier`).
- `surface_depth_in_mm` (float, default `3`) — depositions closer than this
  distance to the n+ surface are stored in the `surface` subtable. It must
  exceed the largest FCCD applied later: deeper depositions count as fully
  active.
- `buffer_len` (str) — LH5 read chunk size (e.g. `"500*MB"`). Controls memory
  usage during processing; does not affect the output.
- `eresmod_default` — energy resolution model applied to non-ON detectors. See
  {ref}`build-tier-hit-hpge` for when this fallback is triggered.
- `aoeresmod_default` — A/E resolution model applied to detectors without a
  per-detector entry. See {ref}`build-tier-hit-hpge` for when this fallback is
  triggered.
- `simulate_psd` (bool, default `True`): enable the single-template A/E PSD
  simulation, written to the `psd/single_temp` subtable of the `hit` tier.
- `simulate_psd_with_psl` (bool, default `False`): enable the
  pulse-shape-library based PSD simulation (see {ref}`hpge-psl-overview`),
  written to the `psd/pulse_lib` subtable. The two flags are independent: enable
  either, both, or neither. Setting both to `False` disables the HPGe PSD
  simulation entirely.
- `tune_hpge_impurities_on_data` (bool, default `False`): select the HPGe
  modeling criteria based on the tuning of the crystal impurity curve to data.
  When `False`, detectors are selected by their operational voltage margin above
  depletion instead. See {ref}`hpge-modeling-criteria`. When `True`, the
  `simlist` items only select the simulations to tune the HPGe impurity model on
  (see {ref}`impurityscan-settings-meta`), and the model file is the only
  default target.
- `psdcuts_default` — PSD cut values applied to detectors without a per-detector
  entry. See {ref}`build-tier-hit-hpge` for when this fallback is triggered.
- `aoemeanmod_default` — A/E mean energy-dependence correction applied to
  detectors without a per-detector entry. Each of the `single_template` and
  `psl` sub-blocks carries an `expression` in the energy `x` (keV) with
  parameters `a`, `b`. The identity model (`a: 0, b: 1`) applies no correction.
  See {ref}`build-tier-hit-hpge` for when this fallback is triggered.
- `two_pass_aoe_correction` (bool, default `False`): when `True`, the `hit` tier
  is built twice. A temporary pre-correction pass computes the raw A/E for the
  simulation IDs matching `simid_regex` (see the `aoemeanmod` par settings, by
  default `sis*_z*_slot*_Pb212_to_Pb208`), from which a per-detector A/E
  energy-dependence correction is fit and applied in the final pass, overriding
  `aoemeanmod_default` where available. The selected experiment (or `simlist`)
  must contain at least one matching `simid`, otherwise the workflow fails at
  DAG-build time.

(evt-tier-settings-meta)=

### `evt` tier

This section specifies how to configure the `evt` tier, which builds physics
events from the hit-level data.

```{code-block} yaml
:caption: simprod/config/tier/evt/{experiment}/settings.yaml

add_random_coincidences: false
random_coincidence_mode: forced_trigger
random_coincidence_runid: null
geds_energy_thr_kev: 25
spms_energy_thr_pe: 0
lar_veto_multiplicity_thr: 4
lar_veto_energy_sum_pe_thr: 4
buffer_len: "50*MB"
skip_opt: false
skip_hit: false
```

- `add_random_coincidences` (bool) — when `true`, random-coincidence (RC) SiPM
  data (taken from `l200data` evt/pet tiers) is mixed in during event building.
- `random_coincidence_mode` (str) — how RC events are selected. With
  `forced_trigger` (the default) the forced/pulser and HPGe-triggered events of
  a physics stream are used. With `noise_trigger` the source is a dedicated
  noise-trigger stream, in which every event is already a random coincidence; no
  trigger flags are read, which is required for configurations without HPGe
  detectors since those carry neither `trigger/is_forced` nor `coincident/*`.
- `random_coincidence_runid` (str | mapping | null) — runid to draw RC data
  from, when it differs from the run being modelled. Defaults to the modelled
  run. Set it when the noise-trigger stream is recorded under its own runid, in
  the format l200-<period>-<run>-<datatype>, e.g. `l200-p13-r003-anp`. It may
  also be a mapping from the modelled runid to its RC runid, for an experiment
  whose runs each have their own noise-trigger stream:

  ```yaml
  random_coincidence_runid:
    l200-p13-r003-aph: l200-p13-r003-anp
    l200-p13-r006-acs: l200-p13-r006-anc
    l200-p13-r008-ath: l200-p13-r008-ant
  ```

  Every modelled run needs its own entry: a missing one raises, rather than
  silently falling back to the run being modelled.

- `geds_energy_thr_kev` (int) — HPGe hit energy threshold in keV; hits below
  this value are discarded.
- `spms_energy_thr_pe` (int) — SiPM hit threshold in photoelectrons; hits below
  this value are discarded.
- `lar_veto_multiplicity_thr` (int) — SiPM multiplicity that marks an event as
  coincident with a LAr signal. See `lar_veto_energy_sum_pe_thr`. LEGEND-200
  uses `4` while LEGEND-1000 uses `1`.
- `lar_veto_energy_sum_pe_thr` (float) — summed SiPM energy in photoelectrons
  that marks an event as coincident with a LAr signal. The `coincident/spms`
  field is `true` when the SiPM multiplicity is `lar_veto_multiplicity_thr` or
  more, or the summed SiPM energy is `lar_veto_energy_sum_pe_thr` or more. Both
  keys are mandatory. LEGEND-200 uses `4` while LEGEND-1000 uses `1`.
- `buffer_len` (str) — LH5 read chunk size (e.g. `"50*MB"`). Controls memory
  usage during processing; does not affect the output.
- `skip_opt` (bool, default `false`) — when `true`, the `opt` (SiPM/LAr) tier is
  skipped: the opt Snakemake rule is not run and the evt output contains only
  HPGe data (no `spms` or `coincident/spms` tables).
- `skip_hit` (bool, default `false`) — when `true`, the `hit` (HPGe) tier is
  skipped: the hit Snakemake rule is not run and the evt output contains only
  SiPM data (no `geds` or `coincident/geds` tables). The hit tier settings and
  `geds_energy_thr_kev` are then optional, and the `par` step builds no HPGe
  parameters: no pulse shape simulation inputs, no modelable-HPGe cache and no
  energy resolution, A/E resolution or PSD cut models.

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

# optional: split the events by the energy deposited in the FCCD layer
fccd_tag_energy_thr_kev: 10

# optional: also build the PDFs on a grid of dead-layer models
morph_dead_layer: true
morph_grid:
  fccd_sigma: [-5, 5, 0.5]
  dlf: [0, 1, 0.1]
morph_min_fccd_in_mm: 0.5
```

- `buffer_len` (str) — LH5 read chunk size (e.g. `"500*MB"`). Controls memory
  usage during processing; does not affect the output.
- `detector_groups` (mapping, optional) — maps group names to Python regex
  strings. Each regex is matched against LEGEND-200 detector names using
  `re.fullmatch`, so `"V.*"` selects all detectors whose names start with `V`.
  When this key is absent, only the implicit `all` group is emitted (equivalent
  to `detector_groups: {all: ".*"}`). Specifying `detector_groups` extends the
  output: every named group is produced in addition to `all`, which is always
  emitted regardless of the config. An explicit `all` entry must be `".*"`; any
  other value is replaced by `".*"` with a warning. See {ref}`pdf-tier` for the
  resulting output schema.
- `fccd_tag_energy_thr_kev` (float, default `null`) — when set, the
  multiplicity-1 events are also split by the energy deposited in the FCCD layer
  of the detector (the depositions closer to the n+ surface than its FCCD).
  Events with more than this energy, in keV, are tagged. The fit can then scale
  them by a separate factor, for example for the rejection of n+ surface events.
  Requires `store_surface_steps` in the `hit` tier (see
  {ref}`hit-tier-settings`).
- `morph_dead_layer` (bool, default `false`) — also build the multiplicity-1
  PDFs on a grid of dead-layer models, in the `pdf_morph` group of the output
  (see {ref}`pdf-tier`). Requires `store_surface_steps` in the `hit` tier.
- `morph_grid` (mapping, default `fccd_sigma: [-5, 5, 0.5]`, `dlf: [0, 1, 0.1]`)
  — the `[first, last, step]` values of the two grid axes. `fccd_sigma` shifts
  the FCCD of each detector by this number of standard deviations of its FCCD
  uncertainty, read from the detector metadata, 0.2 mm if missing. `dlf` is the
  dead layer fraction.
- `morph_min_fccd_in_mm` (float, default `0.5`) — smallest physical FCCD in mm.
  The shifted FCCD is limited to the range from this value to
  `surface_depth_in_mm` of the `hit` tier.

## `pars/` — simulation parameters

Metadata is organized in this directory by experimental configuration (first
level) and detector type (second level), mirroring the `tier/` structure.

The operational voltages of the HPGe detectors are hardware data and are not
stored here. The Simflow reads them from `hardware/configuration/opvs/` in the
metadata (see {func}`legendsimflow.aggregate.gen_hpge_modeling_status`).

(ssd-settings-meta)=

### Pulse shape simulation settings

A single shared YAML file (applies to all detectors and voltages) that sets the
[`SolidStateDetectors.jl`](https://juliaphysics.github.io/SolidStateDetectors.jl/stable/)
simulation control parameters for the
[`build_hpge_drift_time_map`](../api/snakemake_rules.md),
[`build_hpge_pulse_shape_library`](../api/snakemake_rules.md) and
[`build_hpge_psl_scan`](../api/snakemake_rules.md) rules. Missing keys, or a
missing file, keep the defaults shown below. The pulse-shape library scan can
use its own grid size, see {ref}`impurityscan-settings-meta`.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/ssd/settings.yaml (values shown are the defaults)

grid_size_in_mm: 0.5
ssd_refinement_limits: [0.2, 0.1, 0.05, 0.02]
padding: 3
```

| Key                     | Type          | Default                  | Description                                                                                                                                                                                                                         |
| ----------------------- | ------------- | ------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `grid_size_in_mm`       | float         | `0.5`                    | Simulation grid spacing in mm. Execution time scales quadratically with `1/grid_size_in_mm`.                                                                                                                                        |
| `ssd_refinement_limits` | list of float | `[0.2, 0.1, 0.05, 0.02]` | SSD adaptive-mesh refinement thresholds. Each entry drives one refinement pass; smaller values give a more accurate electric field at higher cost. **Overly coarse values can prevent full detector depletion — change with care.** |
| `padding`               | int           | `3`                      | Number of pixel layers padded around the simulated map (drift-time map and pulse-shape library) boundary to avoid grid edge effects.                                                                                                |

:::{tip}

In test or CI environments, setting `grid_size_in_mm: 10.0` reduces the number
of grid points by a factor of ~400 compared to the 0.5 mm production default,
cutting script runtime from many minutes to seconds.

:::

(psl-scan-settings-meta)=

### Pulse shape scan settings

The grid scanned by the pulse shape library scan, which simulates a detector
over a range of impurity profiles and depletion voltages. Each key is a Julia
range written `"start:step:stop"`, quoted so that it is read as a string. When
the file is absent the script falls back to built-in defaults.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/ssd/scan_settings.yaml

depv_shift: "-1000:20:-20"
slope: "-0.9:0.2:3.0"
```

| Key          | Type  | Default          | Description                                                                                                                                                                                                                      |
| ------------ | ----- | ---------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `depv_shift` | range | `"-1000:20:-20"` | Depletion voltages to simulate, in V, relative to the operational voltage of the detector. The impurity density is rescaled to reach each of them.                                                                               |
| `slope`      | range | `"-0.9:0.2:3.0"` | Scaling factors applied to the non-constant part of the impurity profile, dimensionless. A value of 0 leaves the profile as stored in the metadata; the impurity density at the seed end of the boule is preserved in all cases. |

The scan cost is the product of the two range lengths, so the defaults above are
$50 \times 20 = 1000$ detector simulations.

`depv_shift` is relative, but the scan file stores the resulting depletion
voltages as absolute values in V. See
{func}`legendsimflow.psl.validate_ssd_scan_grid` for the layout of that file.

(eresmod-metadata-dir)=

### Energy resolution model defaults

An optional validity-based metadata directory providing HPGe-specific energy
resolution parameters. When present, it can supplement or fully replace
`l200data` as the source of energy resolution parameters — enabling simulations
for experiments that have not yet collected data (e.g. LEGEND-1000). The
structure follows the same validity-based format as
`pars/{experiment}/geds/skip/`.

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

The HPGe A/E mean energy-dependence correction is configured through the
`aoemeanmod_default` key of the {ref}`hit-tier-settings` (applied to every
detector) and, optionally, computed per detector by the Simflow itself when the
`two_pass_aoe_correction` hit-tier setting is enabled.

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
t0_field: auto
end_time_field: geds/psd/low_aoe/time
drift_time_slices: "1000:200:2000"
tier_name: auto
max_files: null
charge_output: wf_pz_win
curr_output: curr_av
energy_output: cuspEmax
```

| Key                 | Type        | Default                 | Description                                                                                                                                                                                                                                                                                                                     |
| ------------------- | ----------- | ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `build_per_runid`   | bool        | `false`                 | When `true`, produce one superpulse LH5 file per `(runid, detector)` instead of one file per detector accumulating all runs. See {ref}`hpge-superpulses-extraction`.                                                                                                                                                            |
| `min_number_wfs`    | int         | `10`                    | Minimum waveforms required to form a superpulse for a slice; slices with fewer events are skipped.                                                                                                                                                                                                                              |
| `target_wfs`        | int         | `100`                   | Target number of waveforms to collect per slice before stopping accumulation.                                                                                                                                                                                                                                                   |
| `chi2_threshold`    | float       | `3`                     | Reduced chi2 threshold for the self-similarity cut; waveforms above this value are rejected before computing the final superpulse.                                                                                                                                                                                              |
| `t0_field`          | str         | `auto`                  | Event-level field used as the start time for drift-time calculation. Events where this field is NaN are discarded; for `spms/event_t0` this removes events without a coincident SiPM signal. With `auto`, use `spms/event_t0` if the files of the run have it, otherwise `geds/t0`; the choice is made separately for each run. |
| `end_time_field`    | str         | `geds/psd/low_aoe/time` | Event-level field used as the end time for drift-time calculation.                                                                                                                                                                                                                                                              |
| `drift_time_slices` | str         | `"1000:200:2000"`       | Drift-time bins as `start:step:stop` in ns; the default creates 200 ns-wide bins from 1000 to 2000 ns.                                                                                                                                                                                                                          |
| `tier_name`         | str         | `auto`                  | Tier to select events from: `evt` or `pet` (event tiers), `hit` or `pht` (hit tiers). With `auto`, use `pet` if the data production has it, otherwise `evt`.                                                                                                                                                                    |
| `max_files`         | int or null | `null`                  | If set, limits the number of raw and evt files processed per run (useful for testing).                                                                                                                                                                                                                                          |
| `charge_output`     | str         | `wf_pz_win`             | DSP processing chain output name for the charge waveform.                                                                                                                                                                                                                                                                       |
| `curr_output`       | str         | `curr_av`               | DSP processing chain output name for the current waveform.                                                                                                                                                                                                                                                                      |
| `energy_output`     | str         | `cuspEmax`              | DSP processing chain output name for the energy estimator used in waveform normalisation.                                                                                                                                                                                                                                       |

`t0_field` and `end_time_field` name fields of the `tier_name` tier: change them
from the defaults when selecting from a `hit` or `pht` tier (`t0_field: auto`
works only with event tiers).

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

(elecmod-fit-settings-meta)=

### Electronics-response fit settings

:::{note}

These settings apply only when the pulse-shape-library (PSL) PSD simulation is
enabled with `simulate_psd_with_psl: True` in {ref}`hit-tier-settings` (default
`False`), and only when the fit actually runs, i.e. when the
{ref}`elecmod-metadata-dir` has no `default` key. When a `default` key is
present the fit is bypassed and these settings are ignored, except `curr_mwa`.
See {ref}`hpge-psl-overview`.

:::

Tuning parameters for the data-driven fit in `extract_electronics_model_pars`
(and in its impurity-curve scan version `extract_elecmod_scan`), which minimises
the mean RMS between simulated and measured current superpulses (across
drift-time slices, via Minuit) to extract the per-detector `sigma` and `tau`.
Settings are read via `get_par_settings(config, "elecmod")` from an optional
`settings.yaml` in the {ref}`elecmod-metadata-dir`. Keys set in the file
override the built-in defaults below; missing keys keep their default.
`extract_elecmod_scan` uses different defaults for three keys:
`tau_start: 30.0`, `tau_limits: [2.0, 100.0]` and `max_calls: 100`. It does not
read `plot_window`.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/elecmod/settings.yaml (values shown are the built-in defaults)

angle: "000"
sigma_start: 10.0
tau_start: 50.0
sigma_limits: [1.0, 50.0]
tau_limits: [1.0, 200.0]
comparison_window: [-500.0, 500.0]
plot_window: [-600.0, 600.0]
weight_power: 2.0
max_calls: 1000
dt_range_tuning: [600.0, 3000.0]
max_num_superpulses: 5
minimiser_mode: "simplex"
curr_mwa:
  length: 48
  num_mw: 3
  mw_type: 0
```

| Key                   | Type           | Default           | Description                                                                                                                                                                                                   |
| --------------------- | -------------- | ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `angle`               | str            | `"000"`           | Crystal-axis angle tag selecting which PSL waveform set (`waveform_{angle}_deg`) is used from the ideal pulse-shape library.                                                                                  |
| `sigma_start`         | float          | `10.0`            | Initial value (ns) for the Gaussian digitizer-bandwidth `sigma`.                                                                                                                                              |
| `tau_start`           | float          | `50.0`            | Initial value (ns) for the preamplifier exponential decay `tau`.                                                                                                                                              |
| `sigma_limits`        | [float, float] | `[1.0, 50.0]`     | Hard `(lo, hi)` bounds (ns) for `sigma` during the fit.                                                                                                                                                       |
| `tau_limits`          | [float, float] | `[1.0, 200.0]`    | Hard `(lo, hi)` bounds (ns) for `tau` during the fit.                                                                                                                                                         |
| `comparison_window`   | [float, float] | `[-500.0, 500.0]` | `(t_min, t_max)` window (ns) relative to the current peak over which the simulation/data RMS is computed; `null` uses the full waveform overlap.                                                              |
| `plot_window`         | [float, float] | `[-600.0, 600.0]` | `(t_min, t_max)` window (ns) relative to the current peak shown in the best-fit diagnostic plots; `null` falls back to `comparison_window`. Does not affect the fit.                                          |
| `weight_power`        | float          | `2.0`             | Data-amplitude weight exponent `p` (`w = \|data\|**p`) applied to the squared residuals before the RMS; `0` is the plain equal-weight RMS, larger values bias the fit toward the current peak and its flanks. |
| `max_calls`           | int            | `1000`            | Maximum number of Minuit function evaluations.                                                                                                                                                                |
| `dt_range_tuning`     | [float, float] | `[600.0, 3000.0]` | Drift-time range (ns); only superpulse slices whose centre falls in this range are used in the fit.                                                                                                           |
| `max_num_superpulses` | int            | `5`               | Maximum number of drift-time slices used in the fit.                                                                                                                                                          |
| `minimiser_mode`      | str            | `simplex`         | Mode of minimisation to use `simplex` `migrad` or `both`                                                                                                                                                      |
| `curr_mwa`            | dict           | see below         | Moving-window average applied to the simulated current, see below.                                                                                                                                            |

|

`curr_mwa` sets the moving-window average that smooths the simulated current
waveforms, like the `curr_av` processor of the data DSP chain smooths the
measured current used for A/E. It is used in the electronics-response fits, in
the drift-time tuning and in the realistic pulse-shape libraries, also when the
fit is bypassed. Keys:

- `length`: window length in ns.
- `num_mw`: number of windows applied in sequence.
- `mw_type`: window type, see `dspeed.processors.moving_window_multi`.

If `l200data` is configured, these settings are ignored (with a warning) and the
values of the `curr_av` processor in the ICPC DSP configuration of `l200data`
are used instead ({func}`legendsimflow.utils.lookup_curr_mw_pars`).

(modeling-settings-meta)=

### HPGe modeling settings

An optional single shared YAML file (applies to all detectors) that tunes the
automatic gate deciding which HPGe detectors are eligible for drift-time map and
current-pulse model generation (see {ref}`hpge-modeling-criteria`). When absent,
the built-in defaults are used. It is only read when
`tune_hpge_impurities_on_data` is `False` in {ref}`hit-tier-settings`.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/modeling/settings.yaml

min_voltage_above_depletion_in_V: 100
```

| Key                                | Type | Default | Description                                                                                                                                                    |
| ---------------------------------- | ---- | ------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `min_voltage_above_depletion_in_V` | int  | `100`   | Minimum margin (in V) by which a detector's operational voltage must exceed its depletion voltage to be considered modelable. Detectors below it are excluded. |

(impurityscan-settings-meta)=

### HPGe impurity-tuning settings

A single shared YAML file (applies to all detectors) that controls the fit of
the HPGe impurity-curve parameters to LEGEND-200 data (`extract_drift_time_scan`
and `extract_hpge_impurity_models`). It is only read when
`tune_hpge_impurities_on_data` is `True` in {ref}`hit-tier-settings`, and the
file must then exist. Missing keys keep the defaults shown below.

```{code-block} yaml
:caption: simprod/config/pars/{experiment}/geds/impurityscan/settings.yaml (values shown are the defaults)

ssd:
  grid_size_in_mm: null
drift_time_scan:
  min_energy_in_keV: 1500
  max_events: null
impurity_fit:
  energy_range_in_keV: [1500, 2500]
  drift_time_uncertainty_in_ns: 50
  waveform_uncertainty: 0.5
  drift_time_observables:
    percentile: 90
    smoothing_in_ns: 50
    min_peak_fraction: 0.25
```

- `ssd.grid_size_in_mm` (float or `null`): spacing in mm of the (r, z) grid of
  the pulse-shape library scan (`build_hpge_psl_scan`). It overrides
  `grid_size_in_mm` of the {ref}`ssd-settings-meta` for this rule only; `null`
  uses that value. The cost of the scan, and of the rules that read it, scales
  quadratically with `1/grid_size_in_mm`.
- `drift_time_scan.min_energy_in_keV` (float): only simulated events depositing
  more than this energy in the HPGe active volume enter the drift-time scan. It
  must not exceed the lower edge of `impurity_fit.energy_range_in_keV`.
- `drift_time_scan.max_events` (int or `null`): maximum number of simulated
  events used per detector and simulation ID; `null` uses all events.
- `impurity_fit.energy_range_in_keV` (`[low, high]`): energy window of the
  events whose drift times are compared between data and simulation.
- `impurity_fit.drift_time_uncertainty_in_ns` (float): uncertainty on the
  drift-time observables in the drift-time $\chi^2$.
- `impurity_fit.waveform_uncertainty` (float): uncertainty on the RMS difference
  between simulated and measured current superpulses (the RMS of
  `extract_elecmod_scan`) in the waveform $\chi^2$.
- `impurity_fit.drift_time_observables`: the two observables of the drift-time
  distribution compared between data and simulation (see
  {func}`legendsimflow.drift_time.drift_time_observables`):
  - `percentile` (float): percentile of the distribution, in percent.
  - `smoothing_in_ns` (float): standard deviation of the Gaussian smoothing
    applied before the peak search.
  - `min_peak_fraction` (float): the first peak is the lowest drift-time peak
    higher than this fraction of the maximum of the smoothed distribution.

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

(predefined-simlists)=

## Predefined simlists

Named lists of simulations are stored in `simprod/config/simlists.yaml`, by
experiment. Items have the `<tier>.<simid>` format of the `simlist`:

```{code-block} yaml
:caption: simprod/config/simlists.yaml

l200cfg09:
  lar_backgrounds:
    - evt.lar_inside
    - evt.lar_outside
```

Select a list in the `simlist` field of the main configuration file with
`~simlists:NAME`, alone or together with other items. On the command line, quote
the value: `--config simlist="~simlists:lar_backgrounds"`.
