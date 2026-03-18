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
- The `.on()` method from `dbetto.TextDB()` involves filesystem queries and can
  be slow. Avoid using it repeatedly in functions invoked at Snakemake DAG build
  time; consider caching strategies instead
- Other conventions are enforced by pre-commit

## Snakemake scripts

Always sanitize the `snakemake` object at the top of every script:

```python
from legendsimflow import nersc

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821
```

- Assign `snakemake` sub-objects to dedicated variables near the start of the
  file for readability
- Global script parameters should use ALL_CAPS variable names
- Large output files should be written to scratch first (`make_on_scratch()`
  from `nersc.py`)
- Add profile logging in disk/compute-intensive scripts (`profile.py`)

## Snakefiles

- Import all Python modules used in the file at the top
- Order rule fields as: `message`, `input`, `params`, `output`, `log`,
  `benchmark`, `script`/`run`/`shell`
- Prefix functions exclusively used in Snakemake rules with `smk_`
- Use `logger.info()` from `snakemake.logging`; never use `print()`
