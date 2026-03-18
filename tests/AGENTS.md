# AGENTS.md — Testing

Tests are stored in `tests/` and managed with Pytest.

- `conftest.py`: fixtures to create mock configuration objects required to test
  package units
- `test_workflow.py`: integration Snakemake testing of the workflow with a dummy
  production (configured in `tests/dummyprod`) that can be tested in CI
- `l200data/`: test data for the LEGEND-200 data production
