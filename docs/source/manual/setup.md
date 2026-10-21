# Installation and configuration

Clone [legend-simflow](https://github.com/legend-exp/legend-simflow) and give it
a custom name:

```console
> git clone git@github.com:legend-exp/legend-simflow <path/to/prod/cycle>
```

We recommend tagging the production cycle with a version number to be used as
folder name (e.g. `path/to/productions/v1.0.0`).

Before a simulation production can be run, the user must configure the run with
a dedicated file and install the required software dependencies (see
{ref}`setup-software`).

## The configuration file

The `simflow-config.yaml` file in the production directory (the root of the
repository) configures the workflow. To start, copy the
[template](https://github.com/legend-exp/legend-simflow/blob/main/templates/default.yaml)
and edit it:

```console
> cp templates/default.yaml simflow-config.yaml
```

Have a look at the other templates in the
[`templates/`](https://github.com/legend-exp/legend-simflow/tree/main/templates)
folder: one of them might be a better starting point for your production.

The defaults of the optional fields, and of the optional fields of the settings
files in the Simflow metadata (see {doc}`meta`), are listed in
{mod}`legendsimflow.defaults`. Any field can be overridden at runtime with
Snakemake's `--config` option (see {ref}`production`).

### What to simulate

`experiment` (required)

: Label of the experiment to simulate. The Simflow metadata uses the same label
for its configuration files. See the
[legend-pygeom-l200](https://legend-pygeom-l200.readthedocs.io) and
[legend-pygeom-l1000](https://legend-pygeom-l1000.readthedocs.io) docs for the
supported labels.

`simlist` (default `all`)

: Simulations to process. Each simulation is labeled by a unique identifier, the
`simid` (e.g. `hpge_bulk_2vbb`), defined in the `simconfig.yaml` files of the
Simflow metadata (see {ref}`simconfig.yaml`). Give them as a list or a
comma-separated string (e.g. `stp.a,hit.b`). `*` or `all` selects all
simulations in the metadata. Each item is `<tier>.<simid>`, a bare `<simid>`
(the last tier in `make_steps` is used) or `~simlists:NAME`, a predefined
simlist (see {ref}`predefined-simlists`).

`simlist_exclude`

: Simulations to remove from `simlist`, same format. A `<tier>.<simid>` item
removes that item only, a bare `<simid>` removes the simulation in every tier.

`runlist`

: LEGEND data taking runs to build pdfs for, as
`<experiment>-<period>-<run>-<type>` (e.g. `l200-p03-r000-phy`) or `~runlists:`
queries. The simulated events are split across these runs by livetime (see
{ref}`run-partitioning`).

`make_steps` (default: all)

: Workflow steps to load, in any order: `vtx`, `stp`, `par`, `opt`, `hit`,
`evt`, `cvt`, `pdf`. Outputs of steps not listed are never rebuilt. A `simlist`
item such as `hit.myid` builds all listed steps up to `hit`. See
{ref}`make-steps` for use cases and caveats.

### Metadata and data production

`metadata_repo` (required)

: `legend-metadata` or `legend1000-metadata`, the repository at
`paths.metadata`. It selects the
{class}`~legendmeta.legendmetadata.LegendMetadata` or
{class}`~legendmeta.legend1000metadata.Legend1000Metadata` class to read the
metadata.

`metadata_version`

: Revision of the metadata repository at `paths.metadata` (see `metadata_repo`)
to check out (anything `git checkout` accepts). Comment it out while developing
metadata.

`validate_l200data` (default `true`)

: Check the data production at `paths.l200data` before any job starts, and abort
if a check fails (see {func}`legendsimflow.validate.l200data_layout` for the
checks). Ignored if `paths.l200data` is not set.

### Paths

Fields of the `paths` block, with the Simflow inputs and outputs. In the
template, `$_` stands for the production directory. The input paths are checked
before any job starts and by `pixi run validate` (see
{func}`legendsimflow.validate.input_paths`).

`generated` (output)

: Root folder of all workflow outputs.

`l200data` (input)

: LEGEND-200 data production cycle (e.g.
`<...>/public/prodenv/prod-blind/ref-v1.0.0`), used to extract parameters such
as the HPGe energy resolution (see {ref}`par-collection-model`) and checked by
`validate_l200data`.

`metadata` (input)

: Clone of [legend-metadata](https://github.com/legend-exp/legend-metadata) or
[legend1000-metadata](https://github.com/legend-exp/legend1000-metadata) (see
`metadata_repo`). Cloned at runtime if missing.

`config` (input)

: Clone of
[legend-simflow-config](https://github.com/legend-exp/legend-simflow-config), a
submodule of legend-metadata. Its content is described in {doc}`meta`.

`optical_maps.lar` (input)

: Liquid argon optical map (photoelectron detection probability map) used in the
`opt` tier (see {ref}`opt-tier-settings`). The Simflow does not produce it.

`optical_maps.lar_patch` (input, optional)

: Map of a smaller region (e.g. around a calibration source) that replaces that
region of `lar`. Either a single map, or a mapping `<simid>: <map>`. Applied by
the `patch_optical_map` rule (see {doc}`../api/snakemake_rules`).

`hpge_impurity_model` (input, optional)

: Tuned HPGe impurity-model file (`pars/hpge/impurity.yaml`) from a separate
impurity-tuning production (see {ref}`impurityscan-settings-meta` and
{ref}`hpge-dtmap-extraction`). Used for the drift-time maps and the ideal
pulse-shape libraries.

`pars` (output)

: Generated parameter files (e.g. parameters extracted from LEGEND-200 data,
geometries, drift-time maps). See {ref}`par-output-files` and
{ref}`par-detinfo`.

`macros` (output)

: Generated _remage_ macros.

`geom` (output, optional)

: Generated geometries (default `{paths.pars}/geom`), built from
{ref}`geom-config.yaml`.

`dtmaps` (output, optional)

: HPGe drift-time maps (default `{paths.pars}/hpge/dtmaps`, see
{ref}`hpge-dtmap-extraction`).

`tier.<tier>` (output)

: Outputs of each tier (e.g. `tier.stp`), named as in {ref}`file-naming`. The
output fields are described in {doc}`output`. Validation plots go to a `plots/`
subfolder (see {ref}`production` to archive them).

`benchmarks` (output)

: Snakemake rule benchmark files.

`log` (output)

: Snakemake rule log files.

### Other fields

`simflow_rng_seed` (default `0`)

: Non-negative integer, combined with the output file name of each _remage_ job
to compute its random seed (the `{SEED}` substitution, see
{ref}`macro-templates-subst`). Jobs are therefore reproducible. Change it to
simulate statistically independent events with the same configuration.

`benchmark.enabled` (default `false`)

: Turn on a benchmarking run (see {ref}`benchmarking`).

`benchmark.n_primaries`

: Primaries per benchmark job, keyed by tier (e.g. `stp: 10_000`).

`runcmd.remage` (default `remage`)

: _remage_ command, e.g. to use a container (see {ref}`setup-software`).

`nersc`

: NERSC I/O options, see {ref}`sites-nersc-io-optim`.

(setup-software)=

## Software dependencies

The first step is obtaining the software, which is fully specified by the
`pyproject.toml` file. Since some dependencies (notably,
[_remage_](https://remage.readthedocs.io)) are only available on
[conda-forge](https://conda-forge.org), we recommend using
[pixi](https://pixi.sh), which can be easily installed following instructions on
the official documentation.

Pixi installs the software the first time you run one of the Simflow pixi tasks
(see {ref}`production`). To load the software in your shell, run `pixi shell`.

<!-- prettier-ignore -->
:::{admonition} For developers only
:class: caution

A plain Python virtual environment, for example with
[uv](https://docs.astral.sh/uv/), can be useful for development:

```console
> uv venv
> source .venv/bin/activate
> uv pip install .
```

It does not provide the non-Python dependencies (e.g. _remage_) nor the pixi
tasks described in {ref}`production`. Do not use it to run a production.

:::

:::{tip}

Do you want to use a _remage_ container instead of the Conda package? Add to the
configuration file, for example:

```yaml
runcmd:
  remage: apptainer run --cleanenv /.../remage_latest.sif
```

:::
