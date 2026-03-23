(production)=

# Production

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
simulations (a "simlist"):

```console
> snakemake --config simlist="mylist.txt" [...]
```

where `mylist.txt` is a text file in the format:

```
stp.fibers_Ra224_to_Pb208
hit.hpge_bulk_2vbb
...
```

One can even just directly pass a comma-separated list:

```console
> snakemake --config simlist="stp.fibers_Ra224_to_Pb208,hit.hpge_bulk_2vbb
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

Find some useful Snakemake command-line options at the bottom of this page.

## Benchmarking runs

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
> snakemake print_benchmark_stats
simid                                runtime [sec]  speed (hot loop) [ev/sec]  evts / 1h  jobs (1h) / 10^8 evts
-----                                -------------  -------------------------  ---------  ---------------------
stp.sis1_z8430_slot2_Bi212_to_Pb208          139.0                     717.70    2583720                     38
stp.sis1_z8580_slot2_Pb214_to_Po214          167.0                     596.99    2149164                     46
stp.sis1_z8630_slot2_Bi212_to_Pb208          135.0                     740.46    2665656                     37
...                                            ...                        ...        ...                    ...
```

Which computes statistics by inspecting the `stp`-tier (_remage_) logs.

:::{note}

The benchmarking statistics refer exclusively to the hot Geant4 simulation loop.
Overheads such as application initialization or remage built-in post processing
are not taken into account.

:::

## Automatically suggesting simulation settings

After completing a benchmarking run, the `make_sim_settings` rule can
automatically compute suggested values for ``primaries_per_job`` and
``number_of_jobs`` that satisfy a set of configurable constraints:

```console
> snakemake make_sim_settings
simid                                                primaries_per_job  number_of_jobs          total_primaries  est. job time [min]
-----                                                -----------------  --------------  ---------------  -------------------
stp.sis1_z8430_slot2_Bi212_to_Pb208                        1,200,000               9           10,800,000                 27.9
stp.sis1_z8580_slot2_Pb214_to_Po214                        1,000,000              10           10,000,000                 27.9
...                                                              ...             ...                 ...                  ...

Constraints: max_runtime=30 min, target_total_primaries=1.00E+07, n_sig_figs=2
Settings written to: generated/sim_settings_suggestion.yaml
```

The constraints are set via the ``settings_opt`` section in the Simflow
configuration file (all fields are optional and have sensible defaults):

```yaml
settings_opt:
  max_runtime_mins: 30        # max job runtime in minutes (default: 30)
  target_total_primaries: 1e7 # minimum total primaries per simid (default: 10_000_000)
  n_sig_figs: 2               # significant figures in the result (default: 2)
```

The rule produces a YAML file (``{paths.generated}/sim_settings_suggestion.yaml``)
that can be directly merged into the ``stp`` tier ``simconfig.yaml``:

```yaml
hpge_bulk_Rn222_to_Po214:
  primaries_per_job: 1200000
  number_of_jobs: 9
sis1_z8430_slot2_Bi212_to_Pb208:
  primaries_per_job: 1200000
  number_of_jobs: 9
...
```

The algorithm:

1. Sets ``primaries_per_job`` to the maximum number of events that fit within
   ``max_runtime_mins`` (rounded *down* to ``n_sig_figs`` significant figures
   to guarantee the time limit is respected).
2. Sets ``number_of_jobs`` to the minimum number of jobs needed so that
   ``primaries_per_job × number_of_jobs ≥ target_total_primaries`` (rounded
   *up* to ``n_sig_figs`` significant figures to guarantee the statistics
   requirement is met).
