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
simid                                             CPU time [ms/ev]  evts / 1h  jobs (1h) / 10^8 evts
-----                                             ----------------  ---------  ---------------------
raw.l200a-birds-nest-K40                                (13s) 2.79    1288475                     77
raw.l200a-birds-nest-Ra224-to-Pb208                   (191s) 38.33      93916                   1064
raw.l200a-fiber-support-copper-Co60                   (223s) 44.69      80558                   1241
...                                                            ...        ...                    ...
```

:::{note}

The CPU time is a good measure of the actual simulation time, since other tasks
(e.g. application loading) are typically not CPU intensive.

:::
