# AGENTS.md — Testing

Python tests are stored in `tests/` and managed with Pytest. Julia tests are in
`workflow/src/LegendSimflow.jl/test/`. Run all tests with `pixi run test`.

- `conftest.py`: fixtures to create mock configuration objects required to test
  package units
- `test_workflow.py`: integration Snakemake testing of the workflow with a dummy
  production (configured in `tests/dummyprod`) that can be tested in CI
- `scripts/`: tests for the tier scripts in
  `workflow/src/legendsimflow/scripts/`, exercising their standalone CLI
  entrypoints
- `l200data/`: test data for the LEGEND-200 data production

## Test data

`tests/dummyprod/inputs/` contains a standalone metadata instance (hardware
detector specs, channelmaps, datasets) committed directly to the repository.

The dummy production uses two experiments:

- `legend`: a generic experiment name used for unit tests and DAG-building
  tests; its runlist uses real p02 run IDs from the metadata but is **not**
  intended to run an actual production
- `l200cfg01`: used by the `test_stp_workflow` integration test, which exercises
  the full vtx→stp pipeline with remage using the public `legend-pygeom-l200`
  geometry; also used by `test_full_workflow` (`needs_nersc` marker), which runs
  the full vtx→cvt pipeline requiring access to `l200data` — run manually at
  NERSC with `pixi run test-full-workflow`

`legend_testdata` (from `legendtestdata`) is still available as a pytest fixture
for tests that require LH5 data files or other binary assets from the testdata
repository (e.g. `test_reboost.py`, `test_hpge_pars.py`).

## Integration tests (`test_workflow.py`)

The three workflow tests form a progression:

1. **`test_dag`** — touch executor, no remage; verifies DAG resolution
2. **`test_stp_workflow`** (`needs_remage`) — runs vtx→stp with real remage
3. **`test_full_workflow`** (`needs_nersc`, `needs_remage`) — full vtx→cvt
   pipeline, requires `l200data`, NERSC-only

Each test uses a separate output directory to avoid Snakemake cache
cross-contamination.

Adding a new test detector requires consistent entries across diodes, crystals,
channelmaps, statuses, and OPV configs. Dummy detector files should be sourced
from `pylegendtestdata` templates (one per type: V, B, C, P), renamed to the
wanted detector name. The real metadata in `./inputs/` can serve as reference
for realistic values and structure.

Snakemake tracks which targets are up-to-date, so there is no need to clean the
generated directory when only a higher tier fails. A full clean is advisable
once in a while to verify the pipeline works end-to-end from scratch.
