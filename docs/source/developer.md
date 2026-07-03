# Developer's guide

This page gives an overview of development practices for the _legend-simflow_
project. It provides guidance on setting up a development environment, running
tests, building the documentation, and contributing code.

## Setting up a development environment

Fork the [legend-simflow](https://github.com/legend-exp/legend-simflow)
repository on GitHub, then follow the instructions in {doc}`/manual/setup` to
clone and install the software.

**pixi** automatically installs the package in editable mode (configured via
`pyproject.toml`), so `pixi shell` already gives you a live development
environment.

When using a Python-only package manager such as
[uv](https://docs.astral.sh/uv/), install in editable mode with all optional
dependencies (tests, docs, pre-commit):

```console
> uv pip install -e ".[all]"
```

## Running the test suite

The test suite consists of Python tests (managed with
[pytest](https://pytest.org)) and Julia tests. Run the full suite with:

```console
> pixi run -e test test
```

To run only the Python tests:

```console
> pixi run -e test test-python
```

To run only the Julia tests:

```console
> pixi run -e test test-julia
```

:::{note}

The tests should be run in the pixi "test" environment to have access to all
test dependencies!

:::

To run a specific Python test or test function:

```console
> pytest tests/test_workflow.py
> pytest tests/test_workflow.py::test_my_function
```

:::{note}

The Python test suite includes both unit tests (in `tests/`) and integration
tests that exercise the Snakemake workflow with a dummy production configured in
`tests/dummyprod`. The dummy production defines three experiments:

- `legend`: a generic experiment used for unit tests and DAG-building tests; its
  runlist contains real p02 run IDs but is not intended to run an actual
  production
- `l1000dsg01`: used by the `test_l1000_workflow` integration test, which
  exercises the vtx→par pipeline; runs in CI without requiring `l200data`.
  Currently uses l200-p03 runs because l1000 hardware metadata is not yet in
  `dummyprod`
- `l200cfg01`: used by the `test_l200_workflow` integration test (`needs_nersc`
  marker), which runs the full vtx→cvt pipeline and requires access to
  `l200data`; run manually at NERSC

Test data for LEGEND-200 is stored in `tests/l200data/`.

:::

Some integration tests require **remage**, which is only available inside the
pixi `test` environment. The table below summarises what each test needs:

| Test                  | Command                                | Needs pixi | Needs NERSC/l200data |
| --------------------- | -------------------------------------- | ---------- | -------------------- |
| DAG tests             | `pytest tests/test_dag.py`             | no         | no                   |
| `test_l1000_workflow` | `pixi run -e test test-l1000-workflow` | **yes**    | no                   |
| `test_l200_workflow`  | `pixi run -e test test-l200-workflow`  | **yes**    | **yes**              |
| unit / script tests   | `pytest tests/`                        | no         | no                   |

**Always run `pixi run -e test test-l1000-workflow` before committing changes
that touch the workflow, scripts, or metadata.** The dry-run DAG tests alone are
not sufficient to catch runtime failures.

## Code style and linting

We use [pre-commit](https://pre-commit.com) to enforce consistent code style and
catch common errors before committing. The hooks are configured in
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

Every commit must pass all pre-commit hooks. CI enforces this requirement on all
pull requests.

:::

## Plotting convention

When creating plots in the workflow scripts, always call
`legendsimflow.plot.decorate(fig)` before writing the figure to disk (for
example, before `fig.savefig(...)` or `pdf.savefig(fig)`).

## Building the documentation

The documentation is built with [Sphinx](https://www.sphinx-doc.org) and written
in [MyST Markdown](https://myst-parser.readthedocs.io). To build it locally:

```console
> cd docs
> make
```

The output is placed in `build/` (that is, `docs/build/` from the repository
root). Open `build/index.html` in a browser to view the result.

The `make` target first auto-generates API reference pages from Python
docstrings and Snakemake rule docstrings (via `make apidoc`), then runs Sphinx.
To only regenerate the API pages without building the full docs:

```console
> make apidoc
```

To remove all generated files and start fresh:

```console
> make clean
```

:::{tip}

The Sphinx build is configured with `-W --keep-going`, which turns warnings into
errors and reports all of them before stopping. Fix all warnings before
submitting a pull request.

:::

### Documentation conventions

- Pages are written in Markdown using [MyST](https://myst-parser.readthedocs.io)
  syntax. Use colon-fence directives for Sphinx admonitions (`:::` + `{note}`,
  `{tip}`, `{warning}`, etc.).
- Use `console` as the language tag for shell commands (enables copy-button and
  consistent rendering).
- Python docstrings follow the
  [NumPy docstring convention](https://numpydoc.readthedocs.io/en/latest/format.html).
- Snakemake rule docstrings should have a one-sentence summary, a blank line,
  and then a longer description and wildcards section.
- Every new page must be referenced in a `toctree` directive to avoid orphan
  page warnings.

## Resource constraints

Jobs run by Snakemake should not use more than **2 GB of memory** each. The
workflow is designed to run many jobs in parallel on a single node, so
individual jobs must stay within this budget. If a script needs more memory, set
explicit `mem_mb` resources for its rule in the Snakemake profile (see the
`nersc-compute` profile for an example with `build_hpge_drift_time_map`).

## Adding a new tier script

Tier scripts live in `workflow/src/legendsimflow/scripts/tier/`. Every script
must be runnable both from Snakemake and directly from the command line (see
{doc}`/manual/prod` for the user-facing documentation of this feature).

The required pattern uses the `snakemake-argparse-bridge` library. Here is the
minimal structure for a new script `tier/foo.py`:

```python
import argparse

import legenddataflowscripts as ldfs
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, utils
from legendsimflow.scripts import log_script_invocation


@snakemake_compatible(
    mapping={
        "input_files": "input",
        "output_file": "output[0]",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(description="Build the foo tier.")
    parser.add_argument("--input-files", nargs="+", required=True)
    parser.add_argument("--output-file", required=True)
    parser.add_argument("--log-file", default=None)
    parser.add_argument(
        "--simflow-config", "--config", dest="simflow_config", required=True
    )
    args = parser.parse_args()

    config = utils.init_simflow_context(args.simflow_config, workflow=None).config
    log = ldfs.utils.build_log(config.metadata.simprod.config.logging, args.log_file)
    log_script_invocation(log, "tier-foo", parser, args)

    # ... script logic ...


if __name__ == "__main__":
    main()
```

After adding the script, also:

1. **Register the Snakemake rule** in the appropriate `workflow/rules/*.smk`
   file, using `script: "../src/legendsimflow/scripts/tier/foo.py"` and
   providing `input`, `output`, `log`, and `config` fields that match the
   `@snakemake_compatible` mapping.
2. **Add a pixi task** in `pyproject.toml` so users can run the script
   standalone:
   ```toml
   tier-foo = { cmd = "python -m legendsimflow.scripts.tier.foo" }
   ```
3. **Add a test** in `tests/scripts/test_tier_foo.py` following the pattern in
   `tests/scripts/test_tier_cvt.py`.

## Warming Numba caches

Numba `@njit(cache=True)` kernels are compiled to machine code and written to a
shared on-disk cache only on their **first call** with a given type signature.
Importing the module merely defines the kernel; it does not compile it. When
many parallel Snakemake jobs hit an uncached kernel at once, they race to write
the same cache file and segfault.

The `legendsimflow.warmup` module (run by the `warmup` pixi task, a `depends-on`
prerequisite of the production and test tasks) sidesteps this by calling every
such kernel once, serially, before the parallel fan-out, so the jobs only ever
read the cache.

:::{important}

Whenever you add code that calls a new lazily-compiled `@njit(cache=True)`
kernel from a rule that runs in parallel (whether defined here or in a
dependency such as `reboost`), add a call to it in
{func}`legendsimflow.warmup.warm_numba_caches`, using the **exact** argument
dtypes seen at run time (Numba caches per type signature). Verify with
`NUMBA_DEBUG_CACHE=1 pixi run warmup` that the kernel's `.nbc`/`.nbi` files get
written.

:::

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
- **No AI authorship in commit messages**: do not include phrases such as
  "generated by AI", "Co-authored-by: Claude", or similar. The commit message
  should read as if written by the human developer.
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
