# Data names and paths

## File naming conventions

The file naming convention for the `stp`, `opt` and `hit` tiers is:

```
{experiment}-{simid}-job_{jobid}-tier_{tier}.{extension}
```

where each label (in curly brackets `{}`) is alphanumeric (including
underscores: `_`). Do not use dashes (`-`) or other characters.

- `experiment` — name or label for the experimental configuration being
  simulated.

- `simid` — stands for "simulation identifier", i.e. a string to uniquely label
  a simulation. Typically specifies the physical process being simulated and the
  experiment's components involved.

- `jobid` — stands for "(simulation) job identifier". It's a zero-padded integer
  that labels independent jobs across which the simulation is split.

- `tier` — the three-character label of the tier. At the moment the simflow
  supports `vtx`, `stp`, `opt` and `hit` tiers.

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

:::{note}

`simid` keys should preferably adopt "snake case" (lowercase letters, digits,
and underscores only).

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
    from the `vtx` tier (see {ref}`vtx-tier-meta`).

    :::{note}

    If vertices are selected as generator, it means that they include vertex
    position _and_ kinematics. In this situation, the `confinement` key (see
    below) is forbidden.

    :::

- `confinement` — one of:
  - `~defines:NAME` to reference a confinement block in {ref}`confinement.yaml`
  - `~volumes.bulk:PATTERN` to confine to physical volumes matching `PATTERN`
  - `~volumes.surface:PATTERN` to sample on the surface of volumes matching
    `PATTERN`
  - a list of the above strings to combine multiple volume patterns
  - `~vertices:NAME` to used the vertex positions similated by the `vtx` tier
    generator `NAME` (see {ref}`vtx-tier-meta`).
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

### `hit` tier

This section specifies how to configure the post-processing of the
[_remage_](https://remage.readthedocs.io/) simulations from the `stp` tier.

#### Run partitioning

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
table in the `hit` file. A new column named `runid` holding the runid of each
event is appended for convenience.

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

which can be overridden for selected simulations by setting the `runlist` field
for the corresponding `simid` in the `hit`-tier `simconfig.yaml` file:

```{code-block} yaml
:caption: /.../simprod/config/tier/hit/simconfig.yaml

hpge_bulk_Rn222_to_Po214:
  runlist:
    - l200-p03-r003-phy
    - l200-p03-r004-phy
```
