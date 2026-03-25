# Developer's guide

This page gives an overview of development practices for the _legend-simflow_
project. It provides guidance on setting up a development environment, running
tests, building the documentation, and contributing code.

## Setting up a development environment

Start by forking the [legend-simflow](https://github.com/legend-exp/legend-simflow)
repository on GitHub, then clone your fork locally:

```console
> git clone git@github.com:<your-username>/legend-simflow
> cd legend-simflow
```

We recommend using [pixi](https://pixi.sh) as the package manager, as it
handles both Python and conda-forge dependencies (including
[_remage_](https://remage.readthedocs.io)) transparently. To set up and enter
the development environment:

```console
> pixi shell
```

Alternatively, you can use a Python-only package manager like
[uv](https://docs.astral.sh/uv/):

```console
> uv venv
> source .venv/bin/activate
> uv pip install -e ".[all]"
```

The `-e` flag installs the package in editable mode so that changes to the
source are immediately effective without reinstalling.

## Running the test suite

The test suite consists of Python tests (managed with
[pytest](https://pytest.org)) and Julia tests. Run the full suite with:

```console
> pixi run test
```

To run only the Python tests:

```console
> pixi run test-python
```

To run only the Julia tests:

```console
> pixi run test-julia
```

To run a specific Python test or test function:

```console
> pytest tests/test_workflow.py
> pytest tests/test_workflow.py::test_my_function
```

:::{note}

The Python test suite includes both unit tests (in `tests/`) and integration
tests that exercise the Snakemake workflow with a dummy production configured
in `tests/dummyprod`. Test data for LEGEND-200 is stored in `tests/l200data/`.

:::

## Code style and linting

We use [pre-commit](https://pre-commit.com) to enforce consistent code style
and catch common errors before committing. The hooks are configured in
`.pre-commit-config.yaml` and include:

- [Ruff](https://docs.astral.sh/ruff/) — fast Python linter and formatter
- [snakefmt](https://github.com/snakemake/snakefmt) — Snakemake file formatter
- [prettier](https://prettier.io) — Markdown, YAML, and JSON formatter
- [codespell](https://github.com/codespell-project/codespell) — spell checker
- [shellcheck](https://www.shellcheck.net) — shell script validator
- [JuliaFormatter](https://github.com/domluna/JuliaFormatter.jl) — Julia code
  formatter

Install the hooks once after cloning:

```console
> pip install pre-commit
> pre-commit install
```

The hooks then run automatically on every `git commit`. To run them manually on
all files:

```console
> pre-commit run --all-files
```

:::{important}

Every commit must pass all pre-commit hooks. CI enforces this requirement on
all pull requests.

:::

## Building the documentation

The documentation is built with [Sphinx](https://www.sphinx-doc.org) and
written in [MyST Markdown](https://myst-parser.readthedocs.io). To build it
locally:

```console
> cd docs
> make
```

The output is placed in `docs/build/`. Open `docs/build/index.html` in a
browser to view the result.

The `make` target first auto-generates API reference pages from Python
docstrings and Snakemake rule docstrings (via `make apidoc`), then runs
Sphinx. To only regenerate the API pages without building the full docs:

```console
> make apidoc
```

To remove all generated files and start fresh:

```console
> make clean
```

:::{tip}

The Sphinx build is configured with `-W --keep-going`, which turns warnings
into errors and reports all of them before stopping. Fix all warnings before
submitting a pull request.

:::

### Documentation conventions

- Pages are written in Markdown using
  [MyST](https://myst-parser.readthedocs.io) syntax. Use colon-fence (`:::`)
  for Sphinx admonitions (e.g. `:::{note}`, `:::{tip}`, `:::{warning}`).
- Use `console` as the language tag for shell commands (enables copy-button and
  consistent rendering).
- Python docstrings follow the
  [NumPy docstring convention](https://numpydoc.readthedocs.io/en/latest/format.html).
- Snakemake rule docstrings should have a one-sentence summary, a blank line,
  and then a longer description and wildcards section.
- Every new page must be referenced in a `toctree` directive to avoid orphan
  page warnings.

## Contributing

### Git workflow

- Fork the repository and open a pull request from a feature branch.
- Write commit messages following the
  [Conventional Commits](https://www.conventionalcommits.org) specification,
  e.g.:
  - `feat: add support for new tier`
  - `fix: correct macro generation for surface simulations`
  - `docs: update developer guide`
  - `refactor: simplify pattern matching logic`
- Every commit must pass linting (`pre-commit run --all-files`) and the test
  suite (`pixi run test`).

### Pull request checklist

Before opening a pull request:

- [ ] The test suite passes: `pixi run test`
- [ ] All pre-commit hooks pass: `pre-commit run --all-files`
- [ ] The documentation builds without warnings: `cd docs && make`
- [ ] New features or bug fixes are covered by tests
- [ ] New public API symbols include NumPy-style docstrings
- [ ] The PR description references any related issues or pull requests
