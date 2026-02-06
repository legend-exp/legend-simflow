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
  the metadata to label the corresponding configuration files.
- `simlist`: list of simulation identifiers to be processed by Snakemake. Can be
  a list of strings or a path to a text file. If `*` or `all`, will process all
  simulations defined in the metadata.
- `runlist`: list of LEGEND data taking runs to build pdfs for, in the standard
  format `<experiment>-<period>-<run>-<type>` (e.g. `l200-p03-r000-phy`)
- `make_tiers`: list the tiers you would like to populate here. This option is
  useful to speed up the DAG generation and avoid accidentally messing up with
  other tiers.
- `legend_metadata_version`: optionally specify a revision (anything that
  `git checkout` accepts) for the _legend-metadata_ instance used by the
  simflow. If you are _developing_ metadata, comment this option.
- `benchmark`: section used to configure a benchmarking run:
  - `enabled`: boolean flag to enable/disable the feature
  - `n_primaries`: number of primary events to be simulated in the lower tiers
    `ver` and `raw`.
- `paths`: customize paths to input or output files.
  - `l200data`: path to a LEGEND-200 data production folder (e.g.
    `<...>/public/prodenv/prod-blind/ref-v1.0.0`) used to extract production
    parameters (e.g. energy resolution)

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

To load the software, just execute `pixi shell`. Alternatively, you may run
Snakemake via `pixi run snakemake ...`.

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
