# AGENTS.md — Code Conventions

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
- Other conventions are enforced by pre-commit

## Snakemake scripts

Every script under `workflow/src/legendsimflow/scripts/` must be runnable both
from Snakemake and directly from the command line. Follow the pattern in
`scripts/tier/cvt.py` and the full checklist in `docs/source/developer.md`.

- Must accept `--simflow-config` (alias `--config`) and `--log-file` (optional)
- Call `log_script_invocation` right after setting up logging
- Large output files: write to scratch first via `make_on_scratch()`
  (`nersc.py`)
- Add a `tier-<name>` pixi task and a test in
  `tests/scripts/test_tier_<name>.py`
- Add profile logging in disk/compute-intensive scripts (`profile.py`)

## Snakefiles

- Import all Python modules used in the file at the top
- Order rule fields as: `message`, `input`, `params`, `output`, `log`,
  `benchmark`, `script`/`run`/`shell`
- Prefix functions exclusively used in Snakemake rules with `smk_`
- Use `logger.info()` from `snakemake.logging`; never use `print()`
