# AGENTS.md

## Project Overview

This repository hosts the workflow to run the Monte Carlo simulation pipeline
for the LEGEND experiment (often called "the Simflow"). Snakemake is used to
orchestrate routines based on the remage simulation ecosystem.

You can learn about the project from the markdown documentation source in
`docs/source/`:

- `index.md`: overview, key concepts
- `manual/meta.md`: configuration metadata, file naming conventions,
  partitioning
- `manual/setup.md`: install the software and main configuration file
- `manual/prod.md`: how to run the workflow
- `manual/sites.md`: computing site specific docs
- `manual/tips.md`: tips and tricks

## Dependencies & Tooling

Python package. See `pyproject.toml` for dependencies and build config,
including Pixi configuration.

## Common Commands

- Install (dev): `uv pip install -e ".[all]"`
- Test: `pytest -vvv`
- Lint/format: `pre-commit run --all-files`
- Build docs: `cd docs && make`

## Architecture

- `workflow`: project source code
  - `src/legendsimflow`: Python package
    - `patterns.py`: filesystem paths used in Snakemake rules
    - `aggregate.py`: functions to build list of files for the Simflow
    - `hpge_pars.py`/`spms_pars.py`: routines to compute the HPGe/SiPM
      parameters to be used in the `hit` and `opt` tiers, based on information
      from the LEGEND-200 data
    - `scripts`: scripts directly used in Snakemake rules
      - `tier/`: scripts used to build the various tiers
      - `plots/`: scripts used to generate validation plots
      - `libjl/`: small Julia library with code shared by the main Julia scripts
      - `init-julia-env.jl`: script used to setup the Julia environment (see
        `Project.toml` in the same folder) used in the Simflow
  - `rules`: Snakemake modules imported by the main `Snakefile`. Rules belonging
    to each tier are organized in separate modules. The `aux.smk` module is used
    for other auxiliary rules.
  - `profiles`: Snakemake workflow profiles with common Snakemake CLI options
  - `Snakefile`: the main Snakefile
- `templates`: template Simflow configuration files

## Testing

- Stored in `tests/` and managed with Pytest.
- `conftest.py` stores fixtures to create mock configuration objects required to
  test the package units
- `test_workflow.py` is dedicated to integrated Snakemake testing od the
  workflow with a dummy production (configured in `tests/dummyprod`) that can be
  tested in CI
- `l200data` stores test data for the LEGEND-200 data production

## Linting

- Formatting/linting is enforced by pre-commit, hooks are listed in
  `.pre-commit-config.yaml`, while further configuration is in `pyproject.toml`

## Documentation

- Sphinx-based, configuration and extensions stored in `conf.py`
- Build by running `make` from within the `docs/` folder. See `Makefile` for
  other useful Make targets

### Conventions

- Source is stored in `docs/source`, user manual further below `manual/`
- A special script (in `tools/`) is used to generate a Snakemake rule reference
  page from the rule docstrings.
- Docs pages are written in Markdown
- Prefer using colon fences (`:::`) instead of 3-backtick blocks for sphinx
  directives like admonitions. Always add an empty line after the start and
  before the end of a `:::` block, otherwise the `prettier` pre-commit hook will
  re-wrap the paragraph and break the rendering
- Always enable code block highlighting when appropriate. Use the `console` code
  label for shell commands (shell prompt delimiter `>`) and outputs
- Document each function and Snakemake rule with in the respective docstring
- Higher level documentation should be added to the user manual

### Python docstring conventions

- Use NumPy-style docstrings
- Do not document the input and output types in the docstring, they are taken
  automatically from the type annotations
- Wrap prose around the same line length as Ruff's default
- Do not wrap the summary line, it should always fit in one line
- Always reference other Python methods with Sphinx Python-domain
  cross-reference roles using the RST syntax
- Cross-reference methods from other packages when appropriate. A list of
  Interphinx mappings is defined in `docs/source/conf.py`
- All new pages must appear in a toctree - orphan pages cause build warnings
- API reference is auto-generated - do not write it by hand

### Snakemake rule docstring conventions

- Start with a one-sentence imperative-style summary of what the rule does, then
  an empty line
- Summarize what the rule script/command is doing
- At the end, say which wildcards are used by the rule
- Use Markdown formatting syntax, as the myst-parser Sphinx extension is used to
  render docstrings (see `docs/source/conf.py` for configuration)

## Code Conventions

- Make sure every Python source file has a short license statement at the top,
  see other existing files
- Try to adhere to the standard scientific
  [Python code conventions](https://learn.scientific-python.org/development)
- Always add type annotations for input arguments and outputs
- Use as generic a type as possible for arguments, and return as specific a type
  as possible. For example, use `Iterable`/`Sequence` or `Mapping` (from
  `collections.abc`) instead of `list` or `dict` for function arguments, if
  appropriate.
- Prefer `dbetto.AttrsDict` over plain `dict` for non-trivial dictionaries that
  are frequently queried, this allows for attribute-style access rather than
  classic key-style access for elements and makes the code more readable
- The `.on()` method from `dbetto.TextDB()` involves filesystem queries and can
  therefore be slow. Avoid using it repeatedly in functions that are invoked at
  Snakemake DAG build time and consider caching strategies
- Other code conventions are enforced by pre-commit (see linting section)

### Snakemake scripts

- Always sanitize the `snakemake` object as in

  ```python
  from legendsimflow import nersc

  args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821
  ```

  to make the inputs read only

- The various objects carried by the `snakemake` object should be assigned to
  dedicated variables near the start of the file for code readability
- Global script parameters should be given ALL_CAPS variable names and
  initialized near the start of the file
- If the script is writing large output files, they should be written on a
  scratch disk first, see `make_on_scratch()` from the `nersc.py` module
- In disk or compute intensive scripts, it is often useful to add profile
  logging statements during executing and at the end of the script (see
  `profile.py` module)

### Snakefiles

- Always import python modules used in the current file at the top
- Order rule fields as in `message`, `input`, `params`, `output`, `log`,
  `benchmark`, `script`/`run`/`shell`
- Prefix names of functions that are exclusively supposed to be used in
  Snakemake rules with `smk_`
- Use `logger` `from snakemake.logging` to print messagges (e.g.
  `logger.info("hello world")`), never use python `print()`.

## Git Workflow

- Commit style: Conventional Commits (https://www.conventionalcommits.org)
- Run pre-commit hooks before committing, every commit must pass linting and
  testing
- PRs require passing CI
- Make sure to include a concise description of the changes in a PR and link the
  relevant related PRs or issues

## Boundaries

- ALWAYS: make sure linting, testing and docs building succeed
