# Production sites

## NERSC

On NERSC you can use Conda/Mamba and Pixi, but you have to make sure that
software is installed on a performant filesystem like `$SCRATCH` or in
`/global/common/software/...`.

:::{warning}

Pixi cannot install software to `/global/common/software`, because this
filesystem is **read-only on compute nodes** and Pixi needs to modify its
environment at runtime. Use `$SCRATCH` for Pixi environments instead (see the
`detached-environments` setting below).

:::

In case you are running the production on CFS (slower filesystem), these
environment variables are a good starting point (customize paths if needed):

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

Then, you can use Conda/Mamba to install Pixi:

```console
> module load python
> mamba create -n my_env -- 'python<3.14' pixi
> conda config --add channels conda-forge
> conda config --set channel_priority strict
```

The last important step is to tell pixi to not store environments in the current
directory (the simflow directory). At the time of writing, this is only possible
by adding this to `.config/pixi/config.toml`:

```toml
detached-environments = "/pscratch/sd/l/.../software/pixi/envs"
```

:::{note}

Environment variable expansion is not supported in the Pixi config at the time
of writing these instructions, so you have to type in the value of `$SCRATCH` by
hand.

:::

Now you can proceed with setting up and running the production workflow, with
e.g. on a compute node:

```console
> pixi run prod --profile nersc-compute
```

Using the provided `nersc-*` profiles is recommended (have a look at them!).

(sites-nersc-io-optim)=

### I/O optimization

On NERSC, the
[Community File System (CFS)](https://docs.nersc.gov/filesystems/community/) I/O
performance can be a bottleneck for I/O-intensive parts of the workflow (for
example, the `stp` tier production). NERSC provides two relevant mitigations.
The first is the
[scratch filesystem](https://docs.nersc.gov/filesystems/perlmutter-scratch/),
based on solid-state disks and offering very high performance. The second is a
faster, read-only mount of `/global/cfs` at
[`/dvs_ro/cfs`](https://docs.nersc.gov/performance/io/dvs/).

For temporary productions, it is recommended to run the full Simflow entirely on
scratch and, if needed, move the data to CFS at the end. Alternatively, if the
simflow is hosted on CFS, it can automatically read input files from DVS and
temporarily write the output of some I/O-intensive jobs to scratch, then move it
to the expected CFS location at completion. These features can be enabled via
the following block in `simflow-config.yaml`:

```yaml
nersc:
  dvs_ro: true
  scratch: $SCRATCH/<SUBFOLDER>
```

Both features can be disabled by setting the corresponding fields to false.

:::{note}

The NERSC profiles already set the Snakemake option `local-storage-prefix` to
stage intermediate files on scratch. This is a Snakemake-native mechanism and is
**independent** of the Simflow's own `nersc.scratch` setting described above. If
you pass `--default-storage-provider fs` to Snakemake (or add
`default-storage-provider: fs` to a profile), Snakemake will use the
`local-storage-prefix` path as the default storage location for all intermediate
files, copying them to the final output path on completion.
`--default-storage-provider fs` and `nersc.scratch` are **mutually exclusive**
and should not be used together.

:::

:::{warning}

These features are implemented manually for each Snakemake rule, so it could be
that some rules are unaffected by them.

:::

### Multi-node execution

:::{warning}

Multi-node execution via `snakemake-nersc` / `pixi run prodnodes` is **highly
experimental**. Use at your own risk and always verify results afterwards.

:::

As of Snakemake v8.30, support for parallel execution across multiple compute
nodes or interaction with job schedulers (such as Slurm) is not well supported.

:::{note}

An experimental profile to interact with the NERSC batch job system is available
in `nersc-compute-slurm`. Unfortunately, specifying rule resources (which is
required for job submission) seems to slow down the DAG generation step by a
lot.

:::

:::{note}

In principle, one could use the `snakemake-executor-plugin-slurm-jobstep` to
prefix each rule command with a `srun` call, which would make it possible to
parallelize the workflow over several nodes. In practice, NERSC discourages from
starting many `srun` instances for performance reasons. As a result, the only
reliable way to run Snakemake is with one instance on a single compute node.

The `snakemake-nersc` executable, exposed by `legend-simflow` offers a way to
parallelize the workflow in some situations over several nodes. The recommended
invocation is:

```console
> pixi run prodnodes NUMBER_OF_NODES
```

which is equivalent to:

```console
> snakemake-nersc --nodes NUMBER_OF_NODES
```

To submit the workflow as a batch Slurm job instead, use:

```console
> pixi run prodsubmit --time HH:MM:SS [--nodes N]
```

where `--time` is the requested wall time (required) and `--nodes` defaults to
`1`.

The program determines the list of simulations (see the `simlist` in
{ref}`production`) that the user wants to process, partitions it in
`NUMBER_OF_NODES` chunks, and spawns a dedicated Snakemake instance for each,
prefixed by the appropriate `srun` call. This is equivalent to something like:

```sh
srun -N1 -n1 snakemake --workflow-profile workflow/profiles/nersc-compute --config simlist=LIST1 [SNAKEMAKE ARGS] &
srun -N1 -n1 snakemake --workflow-profile workflow/profiles/nersc-compute --config simlist=LIST2 [SNAKEMAKE ARGS] &
...

wait
```

:::{warning}

The Snakemake instances spawned per chunk are fully independent and share no
locking mechanism. It is the user's responsibility to ensure that each instance
operates on a **disjoint subset of the DAG** — otherwise multiple instances may
try to build the same output concurrently, leading to race conditions and
corrupted files.

In practice this means that any shared steps (e.g. `par`, which produces
drift-time maps and other parameters consumed by all `hit`/`opt` jobs) **must be
completed before** launching the multi-node run. The safest approach is to first
run `pixi run prod` with only the shared steps in `make_steps`, wait for it to
finish, and only then launch the multi-node run with those steps excluded.

:::

This approach makes it unfortunately harder to manually interrupt the Simflow,
e.g. hitting `Ctrl+C` will just make Slurm print some jobset status information.
You should instead send signals (`TERM` to stop scheduling more jobs and just
wait for running jobs and `INT` to kill all running jobs) directly to the
snakemake instance.

:::{todo}

Add commands to send signals.

:::

:::{note}

Since the actual jobs that need to be run will not be known a priori, the
_a-priori_ partitioning might be inefficient. To mitigate this, the `simlist` is
randomly shuffled before partitioning.

:::
