# Production sites

## NERSC

On NERSC you can use Conda/Mamba and Pixi, but you have to make sure that
software is installed on a performant filesystem like `$SCRATCH` or in
`/global/common/software/...`. In case you are running the production on CFS
(slower filesystem), these environment variables are a good starting point
(customize paths if needed):

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

Now you can proceed with setting up and running the production workflow.

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

:::{warning

These features are implemented manually for each Snakemake rule, so it could be
that some rules are unaffected by them.

:::
