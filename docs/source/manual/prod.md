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
simulations:

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
> snakemake -q all print_benchmark_stats
simid                                       CPU time [ms/ev]  evts / 1h  jobs (1h) / 10^8 evts
-----                                       ----------------  ---------  ---------------------
stp.birds-nest-K40                                (13s) 2.79    1288475                     77
stp.birds-nest-Ra224-to-Pb208                   (191s) 38.33      93916                   1064
stp.fiber-support-copper-Co60                   (223s) 44.69      80558                   1241
...                                                       ...        ...                    ...
```

:::{note}

The CPU time is a good measure of the actual simulation time, since other tasks
(e.g. application loading) are typically not CPU intensive.

:::
