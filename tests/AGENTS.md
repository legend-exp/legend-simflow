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

The dummy production uses three experiments:

- `legend`: a generic experiment name used for unit tests and DAG-building
  tests; its runlist contains real p02 run IDs but is not intended to run an
  actual production
- `l1000dsg01`: used by `test_l1000_workflow`, which exercises the vtx→pdf
  pipeline; runs in CI without requiring `l200data`. Currently uses l200-p03
  runs (`l200-p03-r000-phy`, `l200-p03-r001-phy`) because l1000 hardware and
  crystal metadata are not yet in `dummyprod`
- `l200cfg01`: used by `test_l200_workflow` (`needs_nersc` marker), which runs
  the full vtx→cvt pipeline requiring access to `l200data` — run manually at
  NERSC with `pixi run -e test test-l200-workflow`

`legend_testdata` (from `legendtestdata`) is still available as a pytest fixture
for tests that require LH5 data files or other binary assets from the testdata
repository (e.g. `test_reboost.py`, `test_hpge_pars.py`).

Large binary files that the Snakemake configs reference (e.g. optical maps) are
**gitignored** and populated at the start of every test session by the
`dummyprod_optmap` autouse fixture in `conftest.py`, which copies them from
`legend_testdata`. Do not commit empty placeholder files for these assets.

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

## Integration tests (`test_workflow.py`)

The remage-driven workflow tests form a progression:

1. **`test_l1000_workflow`** (`needs_remage`) — runs vtx→pdf with real remage,
   experiment `l1000dsg01`; runs in CI. **Requires pixi** (remage is only in the
   pixi environment): `pixi run -e test test-l1000-workflow`. The `l1000dsg01`
   hit settings enable the `simulate_psd_with_psl` tier setting (the live name
   of what used to be the unread `has_detailed_psd` key), so this test also
   exercises the PSL-based "detailed" PSD path: the realistic pulse-shape
   library is built in the par tier, consumed by the hit tier into a `psd_psl`
   sub-table, and read back by the evt tier into `geds/psd_psl`.
2. **`test_l200_workflow`** (`needs_nersc`, `needs_remage`) — full vtx→cvt
   pipeline, experiment `l200cfg01`, requires `l200data`, NERSC-only. Run with:
   `pixi run -e test test-l200-workflow`

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
