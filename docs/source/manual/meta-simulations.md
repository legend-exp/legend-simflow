# Defining the simulations

The simulations of a production are defined per experimental configuration in
`simprod/config/tier/stp/{experiment}/`. Which of them a production runs is set
by the `simlist` of the configuration file (see [](setup.md)).

(simconfig.yaml)=

## `simconfig.yaml`

`tier/stp/{experiment}/simconfig.yaml` maps each simulation identifier (simid)
to the description of the simulation, run by the `stp` step (see [](stp.md)):

```{code-block} yaml
:caption: tier/stp/{experiment}/simconfig.yaml

hpge_bulk_Rn222_to_Po214:
  template: $_/template.mac
  generator: ~defines:Rn222_to_Po214
  confinement: ~defines:hpge_bulk
  primaries_per_job: 10_000
  number_of_jobs: 4
```

A simid may contain letters, digits, underscores and hyphens, but no dots.
Prefer underscores: hyphens separate the fields of the output file names (see
{ref}`file-naming`).

`template` (required)

: the _remage_ macro template of the simulation (see
{ref}`macro-templates-subst`). `$_` stands for the folder of `simconfig.yaml`.

`generator` (optional)

: the primary particles: `~defines:NAME`, a block of {ref}`generators.yaml`, or
`~vertices:NAME`, a vertex generator of the `vtx` step (see
{ref}`vtx-simconfig.yaml`).

`confinement` (optional)

: where the primaries start, one of: `~defines:NAME`, a block of
{ref}`confinement.yaml`; `~volumes.bulk:PATTERN`, uniformly in the physical
volumes whose names match the regular expression `PATTERN`;
`~volumes.surface:PATTERN`, on the surface of those volumes; a list of
`~volumes.bulk:` and `~volumes.surface:` items, combined; `~vertices:NAME`, the
positions of a vertex generator of the `vtx` step;
`~function:module.function(<...>, args)`, a Python function that returns the
macro commands, see
{func}`legendsimflow.commands.get_confinement_from_function`.

`primaries_per_job` (required)

: the number of primaries simulated by each job.

`number_of_jobs` (required)

: the number of independent jobs the simulation is split into.

`macro_substitutions` (optional)

: values for placeholders of the template, as `KEY: value`, used as `$KEY` in
the template.

`geom_config_extra` (optional)

: changes to the geometry configuration for this simulation (see
{ref}`geom-config.yaml`). Each top-level key replaces the one of the geometry
configuration whole: a nested block is not merged.

A vertex generator used as `generator` must give the kinematics of the
primaries. If it gives only the kinematics, the simulation needs a confinement
(not `~vertices:`) for the positions; if it gives the positions too, the
simulation must not have one. A vertex generator used as `confinement` must give
positions, and cannot be combined with a `~vertices:` generator.

YAML anchors and merge keys can be used to share fields between simulations.

(generators.yaml)=

## `generators.yaml`

`tier/stp/{experiment}/generators.yaml` defines the primary particles, as named
lists of _remage_ macro commands, referenced as `~defines:NAME`:

```{code-block} yaml
:caption: tier/stp/{experiment}/generators.yaml

K40:
  - /RMG/Generator/Select GPS
  - /gps/particle ion
  - /gps/energy 0 eV
  - /gps/ion 19 40
```

Each block must be a list of commands.

(confinement.yaml)=

## `confinement.yaml`

`tier/stp/{experiment}/confinement.yaml` defines where the primaries start, as
named lists of _remage_ macro commands, referenced as `~defines:NAME`:

```{code-block} yaml
:caption: tier/stp/{experiment}/confinement.yaml

hpge_bulk:
  - /RMG/Generator/Confine Volume
  - /RMG/Generator/Confinement/Physical/AddVolume V.*
  - /RMG/Generator/Confinement/Physical/AddVolume B.*
```

See the
[_remage_ manual](https://remage.readthedocs.io/en/stable/manual/confinement.html)
for the confinement commands.

(macro-templates-subst)=

## Macro templates

A template is a _remage_ macro with placeholders. The Simflow fills these once
per simulation:

- `$GENERATOR` and `$CONFINEMENT`, with the commands of the `generator` and
  `confinement` fields;
- `$KEY`, with the value of `KEY` in `macro_substitutions`.

A placeholder without a value is an error, so a template that contains
`$CONFINEMENT` needs a simulation with a `confinement`. These are filled for
each job:

- `{N_EVENTS}`, the number of primaries of the job;
- `{SEED}`, the random seed of the job. It is computed from the name of the
  output file of the job, so it is different for every job and the same when a
  job runs again; the `simflow_rng_seed` of the configuration file (see
  [](setup.md)) shifts all the seeds.

```{code-block}
:caption: template.mac

/RMG/Manager/Randomization/Seed {SEED}
...
$GENERATOR
$CONFINEMENT
...
/run/beamOn {N_EVENTS}
```

(geom-config.yaml)=

## Geometry

`geom/{experiment}-geom-config.yaml` configures the geometry generator of the
experiment, which writes the GDML geometry of each simulation (see [](stp.md)).
The `executable` key (required) names the generator, e.g. `legend-pygeom-l200`;
all other keys are options of the generator, see its documentation (e.g.
[legend-pygeom-l200](https://legend-pygeom-l200.readthedocs.io)). `$_` stands
for the folder of the file. A simulation can change the configuration with
`geom_config_extra`.

(vtx-tier-meta)=

(vtx-simconfig.yaml)=

## Vertex generators

`tier/vtx/{experiment}/simconfig.yaml` defines the external vertex generators
used by the `vtx` step (see [](vtx.md)), referenced from a simulation as
`~vertices:NAME`:

```{code-block} yaml
:caption: tier/vtx/{experiment}/simconfig.yaml

hpge_surface:
  command: >-
    revertex hpge-surf-pos --detectors "[VB]*" --surface-type nplus
    --gdml {GDML_FILE} --out-file {OUTPUT_FILE} --n-events {N_EVENTS}
  product: positions
```

`command` (required)

: the shell command that generates the vertices. The Simflow fills `{GDML_FILE}`
(the geometry of the simulation), `{OUTPUT_FILE}` and `{N_EVENTS}` (the
primaries per job of the simulation). Other braces in the command must be
doubled, `{{...}}`, and shell patterns quoted.

`product` (required)

: what the generator writes: `positions`, `kinematics`, or
`kinematics_and_positions`. It decides how a simulation can use the generator,
see above.

(runlist-meta)=

## Run lists

The run list sets the data taking runs among which the simulated events are
partitioned (see {ref}`run-partitioning`). A run is written
`<experiment>-<period>-<run>-<datatype>`, e.g. `l200-p03-r001-phy`. A run list
can also contain queries of the `runlists.yaml` file of
[legend-datasets](https://github.com/legend-exp/legend-datasets), written
`~runlists:<group>.<datatype>.<period>`:

```
~runlists:valid.phy.p04 -> [l200-p04-r000-phy, l200-p04-r001-phy, ...]
```

The runs are sorted and duplicates removed. The run list of the configuration
file (`runlist`, see [](setup.md)) applies to all simulations. A simulation can
have its own, in `tier/hit/{experiment}/simconfig.yaml`:

```{code-block} yaml
:caption: tier/hit/{experiment}/simconfig.yaml

hpge_bulk_Rn222_to_Po214:
  runlist:
    - l200-p03-r003-phy
    - ~runlists:valid.phy.p04
```

It is used by the `opt` and `evt` steps too. A simulation without its own run
list and without a global one is an error.

(predefined-simlists)=

## Simulation lists

`simlists.yaml` names lists of simulations, per experimental configuration:

```{code-block} yaml
:caption: simprod/config/simlists.yaml

l200cfg09:
  lar_backgrounds:
    - evt.lar_inside
    - evt.lar_outside
```

The items have the format of the `simlist` of the configuration file. A list is
selected with `~simlists:NAME` in `simlist` or `simlist_exclude`, alone or
together with other items (see [](setup.md)). On the command line, quote it:
`--config simlist="~simlists:lar_backgrounds"`.
