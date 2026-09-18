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

The metadata comes from two places:

- `tests/dummyprod/inputs/`: a standalone metadata instance (hardware detector
  specs, channelmaps, datasets, simflow configuration) committed directly to
  this repository. It serves the **unit tests** and the tier-script tests.
- the mock `legend-metadata` in
  [legend-testdata](https://github.com/legend-exp/legend-testdata), fetched by
  `legendtestdata`. It holds the hardware and the datasets of the **full-chain
  workflow test**.

The `dummyprod_testdata` fixture assembles what that test reads, in
`tests/dummyprod/legend-metadata`: the hardware and datasets above, next to the
`simprod` configuration committed here, which changes far more often.
`legend-pygeom-l200` builds the geometry from it, always the same one: two
strings of four germanium channels, every one of them the test ICPC `V99999Z`,
surrounded by the fiber shroud and its 58 SiPM channels. Only `V00001A` and
`V00001B` are `on`, so only those two are modelled. The test simulates both the
pulse shape discrimination and the optics.

The dummy production uses two experiments:

- `legend`: a generic experiment name used for unit tests and DAG-building
  tests; its runlist contains real p02 run IDs but is not intended to run an
  actual production
- `l200cfg01`: the LEGEND-200 experiment. `test_l200_workflow` runs it in CI on
  the mock metadata; `test_l200_nersc_workflow` (`needs_nersc`) runs the same
  experiment against `l200data` at NERSC

`legend_testdata` (from `legendtestdata`) is still available as a pytest fixture
for tests that require LH5 data files or other binary assets from the testdata
repository (e.g. `test_reboost.py`, `test_hpge_pars.py`).

Large binary files that the unit-test configs reference (e.g. optical maps) are
**gitignored** and populated at the start of every test session by the
`dummyprod_optmap` autouse fixture in `conftest.py`, which copies them from
`legend_testdata`. Do not commit empty placeholder files for these assets. The
full-chain config reads them straight from the linked testdata instead.

`tests/scripts/conftest.py` is distinct from `tests/conftest.py`. It contains
session-scoped integration fixtures that build the full vtx→cvt pipeline step by
step and cache the outputs for the duration of the test session. The fixtures
are shared across all tests under `tests/scripts/`: `legend_gdml_path`,
`legend_stp_path`, `legend_dtmap_path`, `legend_opt_path`, `legend_hit_path`,
`legend_evt_path`, `legend_cvt_path`.

A pre-built static drift-time map is committed at
`tests/dummyprod/inputs/simprod/V05261B-4200V-hpge-drift-time-map.lh5`. It
contains constant 1000 ns drift times on a 1 mm grid for detector V05261B at
4200 V. The `legend_dtmap_path` fixture uses this file directly so that the
Julia drift-time map script does not need to run during unit tests.

## DAG tests (`test_dag.py`)

Assert on the resolved DAG structure (via `dag.jobs`, not run logs); no
remage/NERSC, run in the default suite. Builds use the **touch** executor on a
throwaway output dir, not a dry run: touch marks the `cache_modelable_hpges`
checkpoint complete, so `smk_load_hpge_cache` falls back to metadata and the
per-detector rules downstream of it (PSL / drift-time map builds) expand (a dry
run leaves them unresolved); the throwaway dir keeps placeholders out of the
real `generated*` dirs.

- `test_dag` / `test_dag_simlist`: full DAG resolves; a simlist target schedules
  the PSD-gated drift-time map plots.
- `test_make_steps_selects_tiers`: `make_steps` selects which tier rules enter
  the DAG (tiers are decoupled, e.g. hit without opt).
- `test_simulate_psd[_with_psl]_toggles_*`: the hit-tier `simulate_psd_with_psl`
  / `simulate_psd` settings (edited in a temp metadata copy) add/remove exactly
  the PSL / drift-time-map rules; guards the YAML-to-DAG wiring the dead
  `has_detailed_psd` key broke.
- `test_skip_{opt,hit}_drops_*` / `..._mutually_exclusive`: the evt-tier
  `skip_opt` / `skip_hit` switches drop the opt / hit jobs (negative case: the
  same `make_steps` is unsatisfiable without the switch); both-skip is rejected
  at build time.

## Full-chain tests (`test_workflow.py`)

- Each test uses a separate output directory and a separate `experiment` to
  avoid Snakemake cache cross-contamination.
- Remember: the `legend` experiment name is reserved for generic unit tests.
- The metadata directory for the experiment to be tested is taken from
  `legend-testdata` (through `pylegendtestdata`), inside `data/metadata`.
- The `legend-metadata/inputs/simprod/config` module is instead stored here, in
  `dummyprod/inputs/simprod/config`. It contains metadata for all tested
  experiments. This part of the metadata stays here because it changes
  frequently, and we want to avoid having to change `legend-testdata` every
  time.
- The `simprod/config` module is injected in the mock test metadata folder taken
  from `legend-testdata`.
- During test execution, Snakemake tracks which targets are up-to-date, so there
  is no need to clean the generated directory when only a higher tier fails. A
  full clean is advisable once in a while to verify the pipeline works
  end-to-end from scratch.

Current full-chain tests:

1. **`test_l200_workflow`** (`needs_remage`) — runs vtx→pdf with real remage on
   the mock `legend-metadata`, experiment `l200cfg01`; runs in CI. **Requires
   pixi** (remage is only in the pixi environment):
   `pixi run -e test test-l200-workflow`.
2. **`test_l200_nersc_workflow`** (`needs_nersc`, `needs_remage`) — full vtx→cvt
   pipeline, same experiment, requires `l200data`, NERSC-only. Run with:
   `pixi run -e test test-l200-nersc-workflow`
