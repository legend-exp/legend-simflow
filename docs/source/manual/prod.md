(production)=

# Production

Run a production by using one of the provided site-specific profiles
(recommended):

```console
> pixi run prod --profile <profile-name>
```

This is equivalent to calling Snakemake directly:

```console
> snakemake --workflow-profile workflow/profiles/<profile-name>
```

If no system-specific profile is needed, omit `--profile` (defaults to
`default`):

```console
> pixi run prod
```

To preview what Snakemake would do without executing anything, use `dry`:

```console
> pixi run dry
```

To mark all output files as up-to-date without rerunning jobs, use `touch`:

```console
> pixi run touch
```

The `--config` command line option is very useful to override configuration
values. It can be used, for example, to restrict the production to a subset of
simulations (a "simlist"):

```console
> snakemake --config simlist="mylist.txt" [...]
```

where `mylist.txt` is a text file in the format:

```text
stp.fibers_Ra224_to_Pb208
hit.hpge_bulk_2vbb
...
```

One can even just directly pass a comma-separated list:

```console
> snakemake --config simlist="stp.fibers_Ra224_to_Pb208,hit.hpge_bulk_2vbb"
```

Remember that Snakemake accepts individual output file paths as arguments. If
supplied, Snakemake will only produce those.

```console
> snakemake /../generated/tier/stp/sis1_z8640_slot3_Pb214_to_Po214/l200p15-sis1_z8640_slot3_Pb214_to_Po214-job_0000-tier_stp.lh5
```

Once the production is over, the `print_stats` rule can be used to display a
table with runtime statistics:

```console
> snakemake -q all print_stats
                                                    wall time [s]         wall time [s]
simid                                                (cumulative)   jobs      (per job)  primaries
-----                                               -------------   ----  -------------  ---------
hit.fibers_Ra224_to_Pb208                                83:20:00    100        0:50:00   1.00E+08
stp.fibers_Ra224_to_Pb208                                58:20:35    100        0:35:00   1.00E+08
stp.fibers_Rn222_to_Po214                                33:20:00    100        0:20:00   1.00E+08
...                                                            ...    ...            ...        ...
```

To archive the validation plots produced during the run into a single tarball,
run the `archive_plots` rule manually after the production is complete:

```console
> snakemake archive_plots
```

This collects all `plots/` subdirectories under `generated/` and packs them into
`tarballs/<cycle>-plots.tar.xz`. The rule is not part of the default `all`
target and must always be triggered explicitly.

Find some useful Snakemake command-line options at the bottom of this page.

## Running tier scripts standalone

Tier scripts can be run directly from the command line, outside of Snakemake.
This is useful to quickly reproduce or debug a failing job: at the start of
every run, each script logs the exact command needed to re-run it manually.

Each supported script has a dedicated pixi task named `tier-<name>`:

```console
> pixi run tier-<name> [options]
```

For example, to run the `cvt` script:

```console
> pixi run tier-cvt --evt-files <evt-files> --cvt-file <cvt-file> \
      --simflow-config <path-to-simflow-config.yaml>
```

Pass `--help` to any task to see its full list of options:

```console
> pixi run tier-cvt --help
```

All scripts share the following conventions:

- `--simflow-config` (alias `--config`) accepts the same configuration file that
  the Snakemake workflow uses as `--configfile`. Variable substitution (e.g.
  `$_`) is performed relative to the directory containing the config file,
  matching the behavior of the Snakemake workflow.
- `--log-file` is optional: if omitted, log output goes to the console.

## Selective step execution with `make_steps`

The `make_steps` configuration option controls which workflow steps are loaded
into the Snakemake DAG. Only the rule files for the listed steps are included,
so outputs of excluded steps are treated as pre-existing source files rather
than targets — Snakemake will neither build nor invalidate them.

**Running a single step in isolation** is the most common use case. For example,
to rebuild only the `hit` tier outputs while leaving `stp` files untouched:

```yaml
make_steps:
  - hit
```

or equivalently at the command line:

```console
> snakemake --config make_steps="[hit]"
```

**Running contiguous steps** works the same way — just list all the steps you
need. To run the full parameter and post-processing chain without re-running the
remage simulations:

```yaml
make_steps:
  - par
  - opt
  - hit
  - evt
  - cvt
```

:::{warning}

The `par` step builds the parameters consumed by both `opt` and `hit` (HPGe
drift-time maps, current-pulse models, energy-resolution and A/E models, and
run-statistics partitioning files). If `opt` or `hit` is included in
`make_steps` without `par`, those parameter files **must already exist on
disk**. Omitting `par` is only safe when the parameters have been produced in a
previous run and have not changed.

:::

:::{note}

`vtx` and `par` cannot be used as step prefixes in `simlist` items (e.g.
`par.myid` is invalid) because neither step produces simid-scoped output files.

:::

## Benchmarking runs

Benchmarking runs measure the simulation speed for each `simid` and produce
suggested `primaries_per_job` and `number_of_jobs` values targeting a ~1-hour
wall time per job and at least 10^7 total events per simid.

### Setup

Create a dedicated production cycle (separate directory) and enable benchmark
mode in the configuration file:

```{code-block} yaml
:caption: simflow-config.yaml

benchmark:
  enabled: true
  n_primaries:
    stp: 10_000   # primaries per benchmark job — enough to get a stable speed estimate

make_steps:
  - stp           # only the remage step is needed for benchmarking
```

The `n_primaries` value controls how many events each benchmark job simulates.
It should be large enough for the remage statistics to stabilize, but small
enough that the job finishes quickly (typically a few minutes). The default
value of 10 000 is a reasonable starting point.

### Running

Run the production as usual. Snakemake will spawn exactly one job per `simid`:

```console
> pixi run prod --profile <profile-name>
```

### Inspecting results

Once all benchmark jobs have completed, summarize the results:

```console
> snakemake -q all print_benchmark_stats
simid                                               runtime [sec]  speed (hot loop) [ev/sec]  evts / 1h  ...rounded  jobs (1h) / 10^8 evts  ...rounded
-----                                               -------------  -------------------------  ---------  ----------  ---------------------  ----------
stp.sis1_z8430_slot2_Bi212_to_Pb208                        139.0                     717.70    2583720     2500000                     38          40
stp.sis1_z8580_slot2_Pb214_to_Po214                        167.0                     596.99    2149164     2000000                     46          50
stp.sis1_z8630_slot2_Bi212_to_Pb208                        135.0                     740.46    2665656     2500000                     37          40
...                                                           ...                        ...        ...         ...                    ...         ...
```

The columns are:

- **`runtime [sec]`** — wall time of the benchmark job as reported by _remage_.
- **`speed (hot loop) [ev/sec]`** — average simulation throughput from the
  _remage_ event-loop statistics (averaged over threads if multithreaded).
- **`evts / 1h`** — exact number of events that would be produced in one hour at
  the measured speed.
- **`...rounded`** — same value rounded down to two significant figures with
  half-integer steps (e.g. 2 583 720 → 2 500 000). This is the suggested
  `primaries_per_job` value.
- **`jobs (1h) / 10^8 evts`** — exact number of 1-hour jobs needed to reach 10^8
  total events.
- **`...rounded`** — same, computed from the rounded `primaries_per_job`.

:::{note}

The speed is extracted from _remage_'s hot Geant4 simulation loop only.
Overheads such as application initialization or _remage_ built-in
post-processing are not included.

:::

### Applying the results

After printing the table, the rule also writes an updated simconfig to
`generated/benchmarks/generated-simconfig.yaml`. This file is a copy of the
source `simconfig.yaml` with `primaries_per_job` and `number_of_jobs` patched
for every simid that has benchmark results:

- `primaries_per_job` is set to the rounded events-per-hour value (suggested ~1
  h per job).
- `number_of_jobs` is computed with ceiling division so that
  `primaries_per_job * number_of_jobs >= 10^7`, and is always at least 1.

YAML anchors, aliases, merge keys, and comments from the original
`simconfig.yaml` are preserved.

To apply the suggested values, inspect the generated file and, if satisfied,
copy it over the source:

```console
> diff generated/benchmarks/generated-simconfig.yaml \
       inputs/simprod/config/tier/stp/<experiment>/simconfig.yaml
> cp  generated/benchmarks/generated-simconfig.yaml \
      inputs/simprod/config/tier/stp/<experiment>/simconfig.yaml
```

:::{note}

Simids defined via YAML merge keys (`<<:`) inherit `primaries_per_job` and
`number_of_jobs` from their anchor. If the benchmarked simid is the
anchor-defining entry, the updated values will propagate to all aliases. If only
the alias entry is benchmarked, an explicit override is added to that entry
only.

:::
