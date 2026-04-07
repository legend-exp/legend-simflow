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

## Test data construction

`conftest.py` checks out a pinned commit of `legend-testdata` and copies the
relevant subdirectories (hardware metadata, datasets, ...) into
`tests/dummyprod/inputs/` using `_copy_skip_existing`, which skips files that
already exist at the destination. Those directories are gitignored, so they
always reflect the checked-out testdata commit and cannot be overridden by
committed files. When updating the testdata commit, verify that the new content
is consistent with the test expectations (channelmap validity dates, detector
lists, dataset entries, ...).
