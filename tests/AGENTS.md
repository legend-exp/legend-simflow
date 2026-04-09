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
- `l200p03`: used by the `test_stp_workflow` integration test, which exercises
  the full vtx→stp pipeline with remage using the public `legend-pygeom-l200`
  geometry; also used by `test_full_workflow` (`needs_nersc` marker), which runs
  the full vtx→cvt pipeline requiring access to `l200data` — run manually at
  NERSC with `pixi run test-full-workflow`

`legend_testdata` (from `legendtestdata`) is still available as a pytest fixture
for tests that require LH5 data files or other binary assets from the testdata
repository (e.g. `test_reboost.py`, `test_hpge_pars.py`).
