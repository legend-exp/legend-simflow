# Installation and configuration

Clone [legend-simflow](https://github.com/legend-exp/legend-simflow) and give it
a custom name:

```console
> git clone git@github.com:legend-exp/legend-simflow <path/to/prod/cycle>
```

We recommend tagging the production cycle with a version number to be used as
folder name (e.g. `path/to/productions/v1.0.0`).

Before a simulation production can be run, the user must configure the run with
a dedicated file and install the required software dependencies.

## The configuration file

The `simflow-config.yaml` file resides in the production directory (the root of
the GitHub repository) and allows to customize the workflow in great detail.
Here's a basic description of its fields:

- `experiment`: labels the experiment to be simulated. The same name is used in
  the Simflow metadata to label the corresponding configuration files. See
  [legend-simflow-config](https://github.com/legend-exp/legend-simflow-config/blob/main/README.md)
  for a list of currently supported experiment labels.
- `simlist`: list of simulation identifiers to be processed by Snakemake. Can be
  a list of strings or a path to a text file. If `*` or `all`, will process all
  simulations defined in the metadata.
- `runlist`: list of LEGEND data taking runs to build pdfs for, in the standard
  format `<experiment>-<period>-<run>-<type>` (e.g. `l200-p03-r000-phy`)
- `make_steps`: list the workflow steps to include in the DAG. Only the
  Snakemake rule files for the listed steps are loaded, so Snakemake will treat
  outputs of excluded steps as pre-existing source files rather than targets to
  (re)build. This is the primary mechanism for protecting existing production
  outputs from accidental re-execution and for speeding up DAG generation.
  Available steps (in canonical order): `vtx`, `stp`, `par`, `opt`, `hit`,
  `evt`, `cvt`. If omitted, all steps are included. The order in the list does
  not matter — steps are always sorted canonically. A few caveats apply:
  - Steps `vtx` and `par` do not produce simid-scoped outputs and cannot appear
    as prefixes in `simlist` items.
  - The `par` step builds the parameters consumed by both `opt` and `hit` (HPGe
    drift-time maps, current-pulse models, energy-resolution models,
    run-statistics partitioning files). If `opt` or `hit` is listed without
    `par`, those parameter files must already exist on disk — otherwise
    Snakemake will fail at runtime trying to locate them.
  - When a `simlist` item such as `hit.myid` is specified, the Simflow builds
    outputs cumulatively for all steps up to `hit` that are present in
    `make_steps`. Steps not listed are skipped even in the cumulative pass.

- `legend_metadata_version`: optionally specify a revision (anything that
  `git checkout` accepts) for the _legend-metadata_ instance used by the
  simflow. If you are _developing_ metadata, comment this option.
- `benchmark`: section used to configure a benchmarking run:
  - `enabled`: boolean flag to enable/disable the feature
  - `n_primaries`: number of primary events to be simulated in the lower tiers
    `ver` and `raw`.
- `paths`: paths to Simflow input or output files (or folder hosting them).
  - `generated` (output): root folder for all workflow outputs
  - `l200data` (input): LEGEND-200 data production cycle (e.g.
    `<...>/public/prodenv/prod-blind/ref-v1.0.0`) used to extract production
    parameters (e.g. energy resolution)
  - `benchmarks` (output): Snakemake rule benchmark files
  - `log` (output): Snakemake rule log files
  - `metadata` (input): Simflow input metadata. This is a clone of the
    [legend-metadata](https://github.com/legend-exp/legend-metadata) repository.
    If not present at runtime, the Simflow will attempt a fresh clone.
  - `config` (input): clone of
    [legend-simflow-config](https://github.com/legend-exp/legend-simflow-config).
    This is distributed as a submodule of legend-metadata.
  - `optical_maps` (dict, input): photoelectron detection probability maps for
    various scintillators used in the `opt` tier. These maps are currently not
    produced by the Simflow and therefore supplied as external input.
    - `lar`: the liquid argon optical map file.
  - `pars` (output): root folder for all generated parameter files (e.g. YAML
    files storing parameters extracted from the LEGEND-200 data, geometry files,
    drift time maps).
  - `macros` (output): generated _remage_ macro files
  - `geom` (optional output): generated simulation geometry files. Defaults to
    `{paths.pars}/geom` if not set.
  - `dtmaps` (optional output): generated HPGe drift time maps. Defaults to
    `{paths.pars}/hpge/dtmaps` if not set.
  - `tier` (dict, output): generated outputs for each tier, keyed by tier name
    (e.g. `tier.stp`, `tier.hit`, ...). Validation plots for each tier are
    stored in a `plots/` subdirectory (e.g. `tier.stp/plots/`).
- `runcmd`: command overrides
  - `remage`: remage command. Useful when not working in a Pixi environment and
    the `remage` command is not available.
- `nersc`: see {ref}`sites-nersc-io-optim`.

:::{tip}

All these configuration parameters can be overridden at runtime through
Snakemake's `--config` option.

:::

For a quick start, just copy over the default configuration file from the
templates:

```console
cp templates/default.yaml simflow-config.yaml
```

and customize it.

## Software dependencies

The first step is obtaining the software, which is fully specified by the
`pyproject.toml` file. Since some dependencies (notably,
[_remage_](https://remage.readthedocs.io)) are only available on
[conda-forge](https://conda-forge.org), we recommend using
[pixi](https://pixi.sh), which can be easily installed following instructions on
the official documentation.

To load the software, just execute `pixi shell`. Alternatively, use the provided
pixi tasks to run the workflow without entering a shell:

- `pixi run prod [--profile <name>]` — run a full production (recommended)
- `pixi run dry` — preview what Snakemake would do without executing
- `pixi run touch` — mark all outputs as up-to-date without re-running jobs
- `pixi run snakemake [ARGS]` — pass arbitrary arguments directly to Snakemake

If you prefer, you can use a Python-only package manager of your choice, for
example [uv](https://docs.astral.sh/uv/):

```console
> uv venv
> source .venv/bin/activate
> uv pip install .
```

:::{tip}

Do you want to use a _remage_ container instead of the Conda package? Add to the
configuration file, for example:

```yaml
runcmd:
  remage: apptainer run --cleanenv /.../remage_latest.sif
```

:::
