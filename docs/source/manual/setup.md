# Installation and configuration

## The configuration file

The `simflow-config.yaml` file in the production directory allows to customize
the workflow in great detail. Here's a basic description of its fields:

- `experiment`: labels the experiment to be simulated. The same name is used in
  the metadata to label the corresponding configuration files.
- `simlist`: list of simulation identifiers (see below) to be processed by
  Snakemake. Can be a list of strings or a path to a text file. If `*` or `all`,
  will process all simulations defined in the metadata.
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

:::{tip}

All these configuration parameters can be overridden at runtime through
Snakemake's `--config` option.

:::

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
