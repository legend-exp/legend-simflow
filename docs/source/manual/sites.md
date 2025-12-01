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
