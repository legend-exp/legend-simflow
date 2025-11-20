# legend-simflow

<img src=".github/logo.jpg" alt="legend-simflow logo" align="left" height="170">

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/legend-exp/legend-simflow?logo=git)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Read the Docs](https://img.shields.io/readthedocs/legend-simflow?logo=readthedocs)](https://legend-simflow.readthedocs.io)
![GitHub issues](https://img.shields.io/github/issues/legend-exp/legend-simflow?logo=github)
![GitHub pull requests](https://img.shields.io/github/issues-pr/legend-exp/legend-simflow?logo=github)
![License](https://img.shields.io/github/license/legend-exp/legend-simflow)

End-to-end Snakemake workflow to run Monte Carlo simulations of signal and
background signatures in the LEGEND experiment and produce probability-density
functions (pdfs). Configuration metadata (e.g. rules for generating simulation
macros or post-processing settings) is stored at
[legend-simflow-config](https://github.com/legend-exp/legend-simflow-config).

## Key concepts

- Simulations are labeled by an unique identifier (e.g. `hpge-bulk-2vbb`), often
  referred as `simid` (simID). The identifiers are defined in
  [legend-simflow-config](https://github.com/legend-exp/legend-simflow-config)
  through `simconfig.yaml` files in tier directories `stp` and `ver`.
- Each simulation is defined by a template macro (also stored as metadata) and
  by a set of rules (in `simconfig.yaml`) needed to generate the actual macros
  (template variable substitutions, number of primaries, number of jobs, etc).
- The production is organized in tiers. The state of a simulation in a certain
  tier is labeled as `<tier>.<simid>`. Snakemake understands this syntax.
- The generated pdfs refer to a user-defined selection of LEGEND data taking
  runs. Such a list of runs is specified through the configuration file.
- The production can be restricted to a subset of simulations by passing a list
  of identifiers to Snakemake.

### Workflow steps (tiers)

1. Tier `ver` building: run simulations that generate Monte Carlo event vertices
   needed to some simulations in the next tier. Simulations that do not need a
   special event vertices will directly start from tier `raw`.
1. Tier `stp` building: run full event simulations. Simulation macro commands
   are generated according to rules defined in the metadata.
1. Tier `hit` building: run the first (hit-oriented) step of simulation
   post-processing. Here, "hit" represents a collection of Geant4 "step" in a
   single detector (sensitive volume). In this tier, hit-wise operations like
   optical map or HPGe detector models application are typically performed. The
   "run partitioning" is also performed at this stage (see below).
1. Tier `evt` building: multiple operations are performed in order to build
   actual events.
1. Tier `pdf` building: summarize `evt`-tier output into histograms (the pdfs).

### Run partitioning

"Run partitioning" refers to incorporating information about the experiment's
data taking runs for which the user wants to build pdfs:

- Partition the `hit` event statistics into fractions corresponding to the
  actual total livetime fraction spanned by each selected run. This information
  is extracted from
  [`legend-metadata/datasets/runinfo.yaml`](https://github.com/legend-exp/legend-datasets/blob/main/runinfo.yaml)
- For each partition, apply HPGe models such as energy resolution or
  pulse-shape.
- ...apply optical models (detection probability lookup tables) for the
  scintillators.
- ...apply detector status flags (available in
  [`legend-metadata/datasets/statuses`](https://github.com/legend-exp/legend-datasets/blob/main/statuses))

## Setting up

### The configuration file

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

> [!TIP]
>
> All these configuration parameters can be overridden at runtime through
> Snakemake's `--config` option.

### Software dependencies

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
> uv pip install '.[runprod]'
```

The `runprod` extra includes dependencies required to actually run the workflow,
like Snakemake.

> [!TIP]
>
> Do you want to use a _remage_ container instead of the Conda package? Add to
> the configuration file, for example:
>
> ```yaml
> runcmd:
>   remage: apptainer run --cleanenv /.../remage_latest.sif
> ```

## Production

Run a production by using one of the provided site-specific profiles
(recommended):

```console
> snakemake --workflow-profile workflow/profiles/<profile-name>
```

If no system-specific profiles are provided, the `--workflow-profile` option can
be omitted. Snakemake will use the `default` profile.

```console
> snakemake
```

The `--config` command line option is very useful to override configuration
values. It can be used, for example, to restrict the production to a subset of
simulations:

```console
> snakemake --config simlist="mylist.txt" [...]
```

where `mylist.txt` is a text file in the format:

```
stp.l200p15-fibers-Ra224-to-Pb208
hit.l200p15-hpge-bulk-2vbb
...
```

One can even just directly pass a comma-separated list:

```console
> snakemake --config simlist="raw.l200a-fibers-Ra224-to-Pb208,hit.l200a-hpge-bulk-2vbb"
```

Remember that Snakemake accepts individual output file paths as arguments. If
supplied, Snakemake will only produce those.

```console
> snakemake /../generated/tier/stp/sis1-z8640-slot3-Pb214-to-Po214/l200p15-sis1-z8640-slot3-Pb214-to-Po214_0000-tier_stp.lh5
```

Once the production is over, the `print_stats` rule can be used to display a
table with runtime statistics:

```console
> snakemake -q all print_stats
                                                          wall time [s]         wall time [s]
simid                                                      (cumulative)   jobs      (per job)  primaries
-----                                                     -------------   ----  -------------  ---------
hit.l200a-fibers-Ra224-to-Pb208                                83:20:00    100        0:50:00   1.00E+08
raw.l200a-fibers-Ra224-to-Pb208                                58:20:35    100        0:35:00   1.00E+08
raw.l200a-fibers-Rn222-to-Po214                                33:20:00    100        0:20:00   1.00E+08
...                                                                 ...    ...            ...        ...
```

Find some useful Snakemake command-line options at the bottom of this page.

### Benchmarking runs

This workflow implements the possibility to run special "benchmarking" runs in
order to evaluate the speed of simulations, for tuning the number of events to
simulate for each simulation run.

1. Create a new, dedicated production cycle (see above)
2. Enable benchmarking in the configuration file and customize further settings
3. Start the production as usual.

Snakemake will spawn a single job (with the number of primary events specified
in the configuration file) for each simulation. Once the production is over, the
results can be summarized via the `print_benchmark_stats` rule:

```console
> snakemake -q all print_benchmark_stats
simid                                             CPU time [ms/ev]  evts / 1h  jobs (1h) / 10^8 evts
-----                                             ----------------  ---------  ---------------------
raw.l200a-birds-nest-K40                                (13s) 2.79    1288475                     77
raw.l200a-birds-nest-Ra224-to-Pb208                   (191s) 38.33      93916                   1064
raw.l200a-fiber-support-copper-Co60                   (223s) 44.69      80558                   1241
...                                                            ...        ...                    ...
```

> [!NOTE]
>
> The CPU time is a good measure of the actual simulation time, since other
> tasks (e.g. application loading) are typically not CPU intensive.

## NERSC-specific instructions

On NERSC, you can use Conda/Mamba to install pixi:

```console
> module load python
> mamba create -n my_env -- 'python<3.14' pixi
> conda config --add channels conda-forge
> conda config --set channel_priority strict
```

But you have to make sure that software is installed on a performant filesystem
like `$SCRATCH` or in `/global/common/software/...`. In case you are running the
production on CFS (slower filesystem), these environment variables are a good
starting point (customize paths if needed):

```bash
export SWPREFIX="$SCRATCH/software"
export XDG_CACHE_HOME="$SCRATCH/cache"

# conda
export CONDA_ENVS_DIRS="$SWPREFIX/conda/envs"
export CONDA_PKGS_DIRS="$XDG_CACHE_HOME/conda/pkgs"

# uv
export UV_PYTHON_INSTALL_DIR="$SWPREFIX/share/uv/python"
export UV_TOOL_DIR="$SWPREFIX/share"
export UV_PYTHON_BIN_DIR="$SWPREFIX/bin"

# NOTE: pixi will still store cache in $XDG_CACHE_HOME
export PIXI_HOME="$SCRATCH/software/pixi"

# numba cache
export NUMBA_CACHE_DIR="$XDG_CACHE_HOME/numba"
```

The last important step is to tell pixi to not store environments in the current
directory (the simflow directory). At the time of writing, this is only possible
by adding this to `.config/pixi/config.toml`:

```toml
detached-environments = "/pscratch/sd/l/.../software/pixi/envs"
```

Note that environment variable expansion is not supported, so you have to type
in the value of `$SCRATCH` by hand.

## Useful Snakemake CLI options

```
usage: snakemake [OPTIONS] -- [TARGET ...]

  --dry-run, --dryrun, -n
                        Do not execute anything, and display what would be done. If you have a very large workflow,
                        use `--dry-run --quiet` to just print a summary of the DAG of jobs. (default: False)
  --workflow-profile WORKFLOW_PROFILE
                        Path (relative to current directory) to workflow specific profile folder to use for
                        configuring Snakemake with parameters specific for this workflow (like resources).
  --jobs, -j N          Use at most N CPU cluster/cloud jobs in parallel. For local execution this is an alias for
                        `--cores` (it is though recommended to use `--cores` in that case). Note: Set to `unlimited`
                        to allow any number of parallel jobs.
  --config, -C [KEY=VALUE ...]
                        Set or overwrite values in the workflow config object. The workflow config object is
                        accessible as variable config inside the workflow. Default values can be set by providing a
                        YAML JSON file (see `--configfile` and Documentation).
  --keep-going, -k      Go on with independent jobs if a job fails during execution. This only applies to runtime
                        failures in job execution, not to errors during workflow parsing or DAG construction.
                        (default: False)
  --rerun-triggers {code,input,mtime,params,software-env} [{code,input,mtime,params,software-env} ...]
                        Define what triggers the rerunning of a job. By default, all triggers are used, which
                        guarantees that results are consistent with the workflow code and configuration. If you
                        rather prefer the traditional way of just considering file modification dates, use `--rerun-
                        trigger mtime`. (default: code input mtime params software-env)
  --rerun-triggers {code,input,mtime,params,software-env} [{code,input,mtime,params,software-env} ...]
                        Define what triggers the rerunning of a job. By default, all triggers are used, which
                        guarantees that results are consistent with the workflow code and configuration. If you
                        rather prefer the traditional way of just considering file modification dates, use `--rerun-
                        trigger mtime`. (default: code input mtime params software-env)
  --force, -f           Force the execution of the selected target or the first rule regardless of already created
                        output. (default: False)
  --forceall, -F        Force the execution of the selected (or the first) rule and all rules it is dependent on
                        regardless of already created output. (default: False)
  --forcerun, -R [TARGET ...]
                        Force the re-execution or creation of the given rules or files. Use this option if you
                        changed a rule and want to have all its output in your workflow updated.
  --rerun-incomplete, --ri
                        Re-run all jobs the output of which is recognized as incomplete. (default: False)
  --list-rules, --list, -l
                        Show available rules in given Snakefile. (default: False)
  --summary, -S         Print a summary of all files created by the workflow. The has the following columns:
                        filename, modification time, rule version, status, plan. Thereby rule version contains the
                        version the file was created with (see the version keyword of rules), and status denotes
                        whether the file is missing, its input files are newer or if version or implementation of
                        the rule changed since file creation. Finally the last column denotes whether the file will
                        be updated or created during the next workflow execution. (default: False)
```
