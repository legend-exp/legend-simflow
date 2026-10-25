# AGENTS.md — Code Conventions

## Naming

- `uid`: detector ID in the simulated geometry (remage/reboost tables, TCM
  `table_key`); `rawid`: channel ID in LEGEND-200 data (`daq.rawid` in channel
  map). Keep distinct even when numerically equal. Exception: output fields
  named after the data format (e.g. evt `geds/rawid`, `spms/rawid`)
- Write "UID" in docs, docstrings, comments and messages

## Python

- License statement at the top of every Python file (see existing files)
- Follow
  [scientific Python conventions](https://learn.scientific-python.org/development)
- Always add type annotations for input arguments and outputs
- Use as generic a type as possible for arguments, and return as specific a type
  as possible. Use `Iterable`/`Sequence` or `Mapping` (from `collections.abc`)
  instead of `list` or `dict` for function arguments, if appropriate.
- Prefer `dbetto.AttrsDict` over plain `dict` for non-trivial dictionaries that
  are frequently queried (enables attribute-style access)
- The `dbetto.TextDB.on()` method involves filesystem queries and can be slow.
  Avoid using it repeatedly in functions invoked at Snakemake DAG build time.
  `LegendMetadata.channelmap()` also calls `.on()`. Consider caching strategies
  instead
- When adding code that needs any lazily-compiled, on-disk-cached Numba kernel
  (`@njit(cache=True)`, a `dspeed`/`reboost` processing chain, or any other
  Numba-cached dependency) from a rule that runs in parallel, warm it in
  `legendsimflow/warmup.py` (`warm_numba_caches`/`warm_hpge_dsp_cache`) with the
  exact runtime dtypes, else parallel jobs race the on-disk cache and segfault.
  Verify with `NUMBA_DEBUG_CACHE=1 pixi run warmup`. See `developer.md`.
- Do not add structure the change does not need:
  - No module constant for a value used once or twice; use a literal or an
    argument default, meaning in the docstring. Exception: user-facing knobs and
    parameters that must stay in sync across callers.
  - No helper function for a few lines used once; inline them.
  - Docstring length follows complexity: one line or none for tiny functions.
    Rationale goes in the commit/PR, not in docstrings or comments.
  - Code comments are typically one or two lines, if needed. They state only
    important-to-remember facts about why this implementation was chosen, that
    can't be deduced by reading the code or cannot be easily verified. Do not
    add comments otherwise.
  - Prefer editing an existing function over adding one.
  - Private functions (name starting with underscore) need only a one-line
    docstring.
  - Tests don't need any docstring most of the time. Add one only if there are
    some facts that are very hard to understand from the test code and are key
    to know.
- Other conventions are enforced by pre-commit
- For plotting code, always call `legendsimflow.plot.decorate(fig)` before any
  `savefig`/`pdf.savefig` call.

## Configuration defaults

- Defaults of optional fields of `simflow-config.yaml` and of the
  `simprod/config` settings files live only in `legendsimflow/defaults.py`
  (`CONFIG`, `TIER_SETTINGS`, `PAR_SETTINGS`). Never write them at the call site
  (no `.get("key", default)`, no `DEFAULT_*` dicts in scripts or rules)
- Read the configuration through `init_simflow_context` (fills `CONFIG`),
  `get_tier_settings`/`get_par_settings` (fill the settings defaults), then
  index directly: `get_tier_settings(config, "hit").simulate_psd`
- New optional field: add it with its default to `defaults.py` and document it
  in `docs/source/manual/` (`tests/test_defaults.py` checks this). A field with
  no meaningful default is required: leave it out
- Optional field whose absence has a meaning: default `None`
- Not covered: `paths` (see `apply_path_defaults`), `simconfig` entries,
  geometry configs, and default values of function arguments

## LEGEND-200 data production (`l200data`)

- Read dataflow metadata (e.g. `dataprod/config`, DSP/hit configs) from
  `<l200data>/inputs/`, never from the Simflow `inputs/` (`config.metadata`):
  the processing must match the data production being analysed. From the Simflow
  `inputs/`, use only `datasets`, `hardware` and `simprod`
- Find tier and `par` paths with `utils.lookup_dataflow_config(l200data)`, not
  hardcoded paths
- Keep `validate.py` in sync with every change to what is read from `l200data`
  (files, folders, LH5 objects and fields, which runs, the settings that switch
  a read on): update the specs (`L200DATA_SPEC`, `DATAFLOW_CONFIG_SPEC`,
  `LH5_SPECS`), `_data_problems` and `tests/test_validate.py`

## Resource constraints

- Jobs run by Snakemake should not use more than **2 GB of memory** each. The
  workflow is designed to run many jobs in parallel on a single node, so
  individual jobs must stay within this budget.

## Snakemake scripts

Every script under `workflow/src/legendsimflow/scripts/` must be runnable both
from Snakemake and directly from the command line. Follow the pattern in
`scripts/tier/cvt.py` and the full checklist in `docs/source/developer.md`.

- Must accept `--simflow-config` (alias `--config`) and `--log-file` (optional)
- Call `log_script_invocation` right after setting up logging
- Decorate `main()` with `@snakemake_compatible` (from
  `snakemake_argparse_bridge`) and supply a `mapping` dict mapping CLI argument
  names to Snakemake object paths (e.g. `"stp_file": "input.stp_file"`); see
  `scripts/tier/cvt.py` for the canonical example
- Large output files: write to scratch first via `make_on_scratch()`
  (`nersc.py`)
- All read-only input file paths must be normalized with
  `nersc.dvs_ro(config, path)` right after the config is loaded. This rewrites
  `/global/...` paths to `/dvs_ro/...` on NERSC when `config.nersc.dvs_ro` is
  enabled, and is a no-op elsewhere. Output paths must **not** be passed through
  `dvs_ro`.
- Add a `tier-<name>` pixi task and a test in
  `tests/scripts/test_tier_<name>.py`
- Add profile logging in disk/compute-intensive scripts
  (`reboost.make_profiler`)

## Snakefiles

- Import all Python modules used in the file at the top
- Order rule fields as: `message`, `input`, `params`, `output`, `log`,
  `benchmark`, `script`/`run`/`shell`
- Prefix functions exclusively used in Snakemake rules with `smk_`
- Each step's `gen_all_tier_<step>` target lists every output of that step read
  by later steps, gated by the same settings as the reading rule. A run with
  only the earlier steps in `make_steps` must leave nothing for later steps to
  rebuild. A `temp()` output is deleted once its consumers in the DAG are done,
  so its consumers must also be in that target (e.g. per-detector `par` files
  and their merged files in `gen_list_of_all_par_outputs`)
- Use `logger.info()` from `snakemake.logging`; never use `print()`

## Detector-info caches (`pars/detinfo/`)

- Modelable HPGes at DAG time: only via `smk_load_hpge_cache()`, never from
  metadata; pass `cache=` to `aggregate` functions
- `params` reading the cache: input `smk_wait_for_hpge_cache`, not the cache
  file (would rerun on any detector change)
- Output depending on cache content beyond its inputs: cache content in `params`
- Metadata not tracked: after updates,
  `--forcerun cache_modelable_hpges cache_detector_usabilities` in a separate,
  earlier call
- Write cache files deterministically (sorted, no `set` order): reruns are
  skipped only on unchanged checksums
