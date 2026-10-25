"""DAG-structure tests for the dummy production.

These tests resolve the Snakemake dependency graph and assert on its structure.
They are built with the *touch* executor: it constructs the graph and creates
empty placeholder outputs without running any real job, and (unlike a dry run)
it marks the modelable-HPGe checkpoint (``cache_modelable_hpges``) complete. That
makes ``smk_load_hpge_cache`` fall back to computing the detector list from
metadata, so the per-detector rules downstream of the checkpoint (the PSL and
drift-time map builds) expand too; a dry run would leave them unresolved.

Every build is pointed at a throwaway output directory, so the placeholder
outputs never touch the real ``generated*`` dirs and the resolved graph is
deterministic regardless of what already exists on disk. No remage, Julia or
NERSC access is required, so they run in the default test suite.
"""

from __future__ import annotations

import contextlib
import io
import os
import shutil
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from functools import partial
from pathlib import Path

import dbetto
import pytest
import yaml
from snakemake import api as smkapi
from snakemake.exceptions import MissingInputException, WorkflowError
from snakemake.executors import local

from legendsimflow import aggregate as agg

dummyprod = Path(__file__).parent / "dummyprod"
default_config = dummyprod / "simflow-config.yaml"
l200_config = dummyprod / "simflow-config-l200.yaml"
l1000_config = dummyprod / "simflow-config-l1000.yaml"
all_cores = os.cpu_count() or 1

# the file-producing rule of each tier, used to assert tier (de)selection
TIER_BUILD_RULE = {
    "vtx": "build_tier_vtx",
    "stp": "build_tier_stp",
    "opt": "build_tier_opt",
    "hit": "build_tier_hit",
    "evt": "build_tier_evt",
    "cvt": "build_tier_cvt",
    "pdf": "build_tier_pdf",
}
TIERS = tuple(TIER_BUILD_RULE)

# rules building the realistic HPGe pulse-shape library for the detailed,
# PSL-based PSD; gated by the hit-tier `simulate_psd_with_psl` setting
PSL_PSD_RULES = {
    "build_hpge_pulse_shape_library",
    "convolve_hpge_ideal_pulse_shape_lib",
    "merge_hpge_realistic_psls",
    "extract_electronics_model_pars",
    "merge_electronics_model_pars",
}

# rules building the drift-time maps and current-pulse model for the standard
# PSD; gated by the hit-tier `simulate_psd` setting
STD_PSD_RULES = {
    "build_hpge_drift_time_map",
    "merge_hpge_drift_time_maps",
    "plot_hpge_drift_time_maps",
    "aggregate_hpge_ssd_modeling_info",
    "extract_current_pulse_model",
    "merge_current_pulse_model_pars",
}

# rules fitting the HPGe impurity profiles to data; gated by the hit-tier
# `tune_hpge_impurities_on_data` setting
IMPURITY_TUNING_RULES = {
    "build_hpge_psl_scan",
    "extract_elecmod_scan",
    "merge_electronics_model_scan_pars",
    "extract_drift_time_scan",
    "merge_hpge_drift_time_scans",
    "extract_hpge_impurity_models",
}

# metadata folder (under dummyprod) of the experiments of the full-chain tests
FULL_CHAIN_METADATA = {
    "l200cfg01": "legend-metadata",
    "l1000dsg01": "legend1000-metadata",
}


def _output_paths(base: Path) -> dict:
    """Redirect every generated output under *base* (a throwaway directory)."""
    gen = base / "generated"
    return {
        "generated": str(gen),
        "benchmarks": str(gen / "benchmarks"),
        "log": str(gen / "log"),
        "pars": str(gen / "pars"),
        "macros": str(gen / "macros"),
        "pdf_releases": str(gen / "releases"),
        "tier": {t: str(gen / "tier" / t) for t in TIERS},
    }


def _metadata_paths(
    base: Path, experiment: str, settings_by_tier: Mapping[str, Mapping]
) -> dict:
    """Copy the dummyprod metadata under *base* with overridden tier settings.

    The `l200cfg01` and `l1000dsg01` experiments run on the metadata assembled
    by the ``dummyprod_testdata`` fixture, the other ones on the metadata
    committed here. ``settings_by_tier`` maps a tier name to the settings keys to
    overwrite for *experiment* (e.g. ``{"evt": {"skip_opt": True}}``). Editing
    the copy keeps the committed metadata untouched.
    """
    source = FULL_CHAIN_METADATA.get(experiment, "inputs")
    meta = base / source
    shutil.copytree(dummyprod / source, meta, symlinks=False)
    for tier, overlay in settings_by_tier.items():
        f = meta / "simprod/config/tier" / tier / experiment / "settings.yaml"
        data = yaml.safe_load(f.read_text())
        data.update(overlay)
        f.write_text(yaml.safe_dump(data))
    return {
        "metadata": str(meta),
        "config": str(meta / "simprod/config"),
        "optical_maps": {"lar": str(meta / "simprod/l200cfg01-optmap-dummy.lh5")},
    }


def overrides(
    base: Path,
    *,
    make_steps: list[str] | None = None,
    simlist: str | None = None,
    experiment: str | None = None,
    settings_by_tier: Mapping[str, Mapping] | None = None,
) -> dict:
    """Build config overrides: a throwaway output dir plus the requested knobs.

    When *settings_by_tier* is given the metadata is copied to *base* and edited,
    so the tier settings take effect without mutating the committed metadata.
    """
    paths = _output_paths(base)
    if settings_by_tier:
        assert experiment is not None
        paths |= _metadata_paths(base, experiment, settings_by_tier)
    cfg: dict = {"paths": paths}
    if make_steps is not None:
        cfg["make_steps"] = make_steps
    if simlist is not None:
        cfg["simlist"] = simlist
    return cfg


def run_snakemake(
    configfile: Path,
    config_overrides: Mapping,
    *,
    workdir: Path = dummyprod,
    executor: str = "touch",
    targets: Iterable[str | Path] = (),
    forceall: bool = False,
    forcerun: Iterable[str] = (),
):
    """Run Snakemake with *executor* and return the resolved DAG."""
    output = smkapi.OutputSettings(verbose=False)
    with smkapi.SnakemakeApi(output) as api:
        wf_api = api.workflow(
            snakefile=dummyprod / "workflow/Snakefile",
            workdir=workdir,
            config_settings=smkapi.ConfigSettings(
                configfiles=(configfile,), config=dict(config_overrides)
            ),
            storage_settings=smkapi.StorageSettings(),
            resource_settings=smkapi.ResourceSettings(cores=all_cores),
        )
        dag = wf_api.dag(
            dag_settings=smkapi.DAGSettings(
                targets={str(t) for t in targets},
                forceall=forceall,
                forcerun=set(forcerun),
            )
        )
        # the touch and dry-run executors dump the full job list to stdout
        with contextlib.redirect_stdout(io.StringIO()):
            dag.execute_workflow(executor=executor)
        # the public API exposes no DAG accessor, so reach into the workflow
        # object for the resolved graph
        return wf_api._workflow.dag


def dag_rule_names(
    configfile: Path, config_overrides: Mapping, targets: set[str] | None = None
) -> set[str]:
    """Return the set of rule names making up the resolved workflow DAG.

    Built with the touch executor (see the module docstring), so rules both
    upstream and downstream of the modelable-HPGe checkpoint appear, independent
    of what exists on disk. *targets* are rule names to build instead of the
    default ``all`` rule. Raises if the DAG cannot be resolved (e.g. an input
    has no producing rule).
    """
    dag = run_snakemake(configfile, config_overrides, targets=targets or ())
    return {job.rule.name for job in dag.jobs}


def test_dag(tmp_path):
    """The default (full) DAG resolves without error and covers every tier."""
    rules = dag_rule_names(default_config, overrides(tmp_path))
    assert "all" in rules
    assert set(TIER_BUILD_RULE.values()) <= rules
    assert "write_l1000_hpge_metadata" not in rules


def test_dag_simlist(tmp_path):
    """An explicit simlist target pulls in the PSD-gated drift-time map plots.

    The par tier contributes the drift-time map plots through
    ``process_simlist``; requesting a hit-tier simid must schedule them.
    """
    rules = dag_rule_names(
        default_config, overrides(tmp_path, simlist="hit.pen_plates_Ra224_to_Pb208")
    )
    assert "build_tier_hit" in rules
    assert "plot_hpge_drift_time_maps" in rules


def test_dag_simlist_reads_hpge_checkpoint(tmp_path, config):
    """A simlist target takes the modelable HPGes from the checkpoint output.

    The per-detector rules read the on-disk checkpoint output, which can differ
    from the metadata (e.g. after a metadata update). Here it marks a detector
    that the metadata calls modelable as not modelable.
    """
    detinfo = agg.pivot_detinfo(agg.gen_list_of_all_hpges_valid_for_modeling(config))
    dets = next(iter(detinfo["is_modelable"].values()))
    dets[next(d for d, flag in dets.items() if flag)] = False

    detinfo_dir = tmp_path / "generated/pars/detinfo"
    detinfo_dir.mkdir(parents=True)
    for flag in ("is_modelable", "operational_voltage_in_V"):
        dbetto.utils.write_dict(detinfo[flag], detinfo_dir / f"{flag}.yaml")

    rules = dag_rule_names(
        default_config,
        overrides(
            tmp_path,
            simlist="hit.pen_plates_Ra224_to_Pb208",
            experiment="legend",
            settings_by_tier={"hit": {"simulate_psd_with_psl": True}},
        ),
    )
    assert "build_tier_hit" in rules


@pytest.mark.parametrize(
    ("make_steps", "present", "absent"),
    [
        # only the lowest tiers
        (["vtx", "stp"], {"vtx", "stp"}, {"opt", "hit", "evt", "cvt", "pdf"}),
        # tiers are decoupled: hit can be built without opt
        (["vtx", "stp", "par", "hit"], {"vtx", "stp", "hit"}, {"opt", "evt"}),
        # the full pipeline
        (
            ["vtx", "stp", "par", "opt", "hit", "evt", "cvt", "pdf"],
            {"vtx", "stp", "opt", "hit", "evt", "cvt", "pdf"},
            set(),
        ),
    ],
)
def test_make_steps_selects_tiers(tmp_path, make_steps, present, absent):
    """`make_steps` controls which tiers (and their build rules) enter the DAG."""
    rules = dag_rule_names(default_config, overrides(tmp_path, make_steps=make_steps))
    assert {TIER_BUILD_RULE[t] for t in present} <= rules
    assert {TIER_BUILD_RULE[t] for t in absent}.isdisjoint(rules)


def test_vtx_only_make_steps_warns(tmp_path, caplog):
    """A vtx-only `make_steps` leaves the default target empty, with a warning."""
    with caplog.at_level("WARNING", logger="snakemake"):
        rules = dag_rule_names(default_config, overrides(tmp_path, make_steps=["vtx"]))

    assert rules == {"all"}
    assert "the default target is empty" in caplog.text


# the PSD switches gate per-detector rules downstream of the modelable-HPGe
# checkpoint; the touch executor (see dag_rule_names) is what expands them
STEPS_TO_HIT = ["vtx", "stp", "par", "opt", "hit"]


@pytest.mark.parametrize(
    ("configfile", "experiment"),
    [(l200_config, "l200cfg01"), (l1000_config, "l1000dsg01")],
)
def test_simulate_psd_with_psl_toggles_psl_rules(tmp_path, configfile, experiment):
    """The hit-tier `simulate_psd_with_psl` setting gates the PSL build chain."""
    on = dag_rule_names(
        configfile,
        overrides(
            tmp_path / "on",
            make_steps=STEPS_TO_HIT,
            experiment=experiment,
            settings_by_tier={"hit": {"simulate_psd_with_psl": True}},
        ),
    )
    off = dag_rule_names(
        configfile,
        overrides(
            tmp_path / "off",
            make_steps=STEPS_TO_HIT,
            experiment=experiment,
            settings_by_tier={"hit": {"simulate_psd_with_psl": False}},
        ),
    )
    assert on >= PSL_PSD_RULES
    assert PSL_PSD_RULES.isdisjoint(off)
    # flipping the switch changes nothing but the PSL rules
    assert on - off == PSL_PSD_RULES


@pytest.mark.parametrize(
    ("configfile", "experiment"),
    [(l200_config, "l200cfg01"), (l1000_config, "l1000dsg01")],
)
def test_par_step_builds_all_hit_pars(tmp_path, configfile, experiment):
    """Adding `hit` to a par-only `make_steps` adds only `hit` tier rules.

    Otherwise a par-only production leaves `par` work for the `hit` step.
    """
    settings = {"hit": {"simulate_psd": True, "simulate_psd_with_psl": True}}
    par = dag_rule_names(
        configfile,
        overrides(
            tmp_path / "par",
            make_steps=["vtx", "stp", "par"],
            experiment=experiment,
            settings_by_tier=settings,
        ),
    )
    hit = dag_rule_names(
        configfile,
        overrides(
            tmp_path / "hit",
            make_steps=["vtx", "stp", "par", "hit"],
            experiment=experiment,
            settings_by_tier=settings,
        ),
    )
    assert par >= PSL_PSD_RULES | STD_PSD_RULES
    assert hit - par == {"build_tier_hit", "plot_tier_hit_observables"}


def test_simulate_psd_toggles_dtmap_rules(tmp_path):
    """The hit-tier `simulate_psd` setting gates the drift-time map build chain."""
    on = dag_rule_names(
        l200_config,
        overrides(
            tmp_path / "on",
            make_steps=STEPS_TO_HIT,
            experiment="l200cfg01",
            settings_by_tier={
                "hit": {"simulate_psd": True, "tune_hpge_impurities_on_data": False}
            },
        ),
    )
    off = dag_rule_names(
        l200_config,
        overrides(
            tmp_path / "off",
            make_steps=STEPS_TO_HIT,
            experiment="l200cfg01",
            settings_by_tier={
                "hit": {"simulate_psd": False, "tune_hpge_impurities_on_data": False}
            },
        ),
    )
    assert on >= STD_PSD_RULES
    assert STD_PSD_RULES.isdisjoint(off)
    # flipping the switch changes nothing but the drift-time map rules
    assert on - off == STD_PSD_RULES


def test_tune_hpge_impurities_builds_tuning_chain(tmp_path):
    """The impurity fit resolves down to the scan simulations and the data.

    The fit output is requested directly: the default target reaches it only
    through the par step, which needs every simid to have a single run. The
    simlist selects the one tuning simid, which has a single run.
    """
    cfg = overrides(
        tmp_path,
        simlist="hit.birds_nest_K40",
        experiment="legend",
        settings_by_tier={
            "hit": {"simulate_psd": True, "tune_hpge_impurities_on_data": True}
        },
    )
    cfg["paths"]["l200data"] = str(Path(__file__).parent / "l200data/v3.0.0")
    rules = dag_rule_names(
        default_config, cfg, targets={"extract_hpge_impurity_models"}
    )
    assert rules >= IMPURITY_TUNING_RULES
    assert {"build_superpulses_from_data", "build_tier_stp"} <= rules


def test_tune_hpge_impurities_default_target(tmp_path):
    """The default target of a tuning production is the impurity model only.

    The simlist items name the tuning simids: their tier is not a request to
    build that tier.
    """
    cfg = overrides(
        tmp_path,
        simlist="pdf.birds_nest_K40",
        experiment="legend",
        settings_by_tier={
            "hit": {"simulate_psd": True, "tune_hpge_impurities_on_data": True}
        },
    )
    cfg["paths"]["l200data"] = str(Path(__file__).parent / "l200data/v3.0.0")
    rules = dag_rule_names(default_config, cfg)
    assert "extract_hpge_impurity_models" in rules
    assert TIER_BUILD_RULE["pdf"] not in rules


def test_skip_opt_drops_opt_tier(tmp_path):
    """The evt-tier `skip_opt` setting drops the opt tier from the DAG.

    Setting ``skip_opt`` removes the opt-tier file from the ``build_tier_evt``
    inputs; with opt left out of ``make_steps`` the opt jobs then have no
    consumer and disappear from the graph. Without ``skip_opt`` the same
    ``make_steps`` is unsatisfiable (evt requires an opt file that no rule
    produces), which is exactly what the switch is there to allow.
    """
    steps = ["vtx", "stp", "par", "hit", "evt"]  # note: no opt
    rules = dag_rule_names(
        l200_config,
        overrides(
            tmp_path / "skip",
            make_steps=steps,
            experiment="l200cfg01",
            settings_by_tier={"evt": {"skip_opt": True}},
        ),
    )
    assert {"build_tier_evt", "build_tier_hit"} <= rules
    assert "build_tier_opt" not in rules

    with pytest.raises(MissingInputException):
        dag_rule_names(
            l200_config,
            overrides(
                tmp_path / "noskip",
                make_steps=steps,
                experiment="l200cfg01",
                settings_by_tier={"evt": {"skip_opt": False}},
            ),
        )


def test_skip_hit_drops_hit_tier(tmp_path):
    """The evt-tier `skip_hit` setting drops the hit tier from the DAG."""
    steps = ["vtx", "stp", "par", "opt", "evt"]  # note: no hit
    rules = dag_rule_names(
        l200_config,
        overrides(
            tmp_path / "skip",
            make_steps=steps,
            experiment="l200cfg01",
            settings_by_tier={"evt": {"skip_hit": True}},
        ),
    )
    assert {"build_tier_evt", "build_tier_opt"} <= rules
    assert "build_tier_hit" not in rules

    with pytest.raises(MissingInputException):
        dag_rule_names(
            l200_config,
            overrides(
                tmp_path / "noskip",
                make_steps=steps,
                experiment="l200cfg01",
                settings_by_tier={"evt": {"skip_hit": False}},
            ),
        )


def test_skip_hit_needs_no_hit_settings(tmp_path):
    """With `skip_hit`, the hit settings and the HPGe threshold may be absent.

    The HPGe pulse-shape rules are dropped, although `simulate_psd` defaults to
    true and par is in `make_steps`.
    """
    cfg = overrides(
        tmp_path,
        make_steps=["vtx", "stp", "par", "opt", "evt"],
        experiment="l200cfg01",
        settings_by_tier={"evt": {"skip_hit": True}},
    )
    tier_dir = Path(cfg["paths"]["config"]) / "tier"
    (tier_dir / "hit/l200cfg01/settings.yaml").unlink()
    evt_settings = tier_dir / "evt/l200cfg01/settings.yaml"
    data = yaml.safe_load(evt_settings.read_text())
    del data["geds_energy_thr_kev"]
    evt_settings.write_text(yaml.safe_dump(data))

    rules = dag_rule_names(l200_config, cfg)
    assert {"build_tier_evt", "build_tier_opt"} <= rules
    assert (STD_PSD_RULES | PSL_PSD_RULES).isdisjoint(rules)


def test_skip_opt_and_hit_are_mutually_exclusive(tmp_path):
    """Skipping both the opt and hit tiers is rejected at DAG-build time."""
    with pytest.raises(WorkflowError, match="skip_opt and skip_hit"):
        dag_rule_names(
            l200_config,
            overrides(
                tmp_path,
                make_steps=["vtx", "stp", "par", "evt"],
                experiment="l200cfg01",
                settings_by_tier={"evt": {"skip_opt": True, "skip_hit": True}},
            ),
        )


def test_geom_plots_scheduled(tmp_path):
    """The geometry validation plots are scheduled with the stp tier."""
    rules = dag_rule_names(default_config, overrides(tmp_path))
    assert "plot_geom_rendering" in rules
    assert "plot_geom_hpge_mass" in rules


def test_l1000_dag(tmp_path):
    """The LEGEND-1000 DAG resolves on the mock `legend1000-metadata`.

    The channel names come from the `legend-pygeom-l1000` configuration. The
    metadata holds only the default records, so the detector records are
    derived from them, and `write_l1000_hpge_metadata` writes the HPGe records.
    Both HPGe detectors are modelable, so the drift-time map and PSL chains are
    scheduled. The HPGe mass plot exists only for `legend-pygeom-l200`.
    """
    rules = dag_rule_names(l1000_config, overrides(tmp_path))
    # the only simulation confines its primaries to volumes, so no vertices
    assert set(TIER_BUILD_RULE.values()) - rules == {"build_tier_vtx"}
    assert rules >= STD_PSD_RULES | PSL_PSD_RULES
    assert "write_l1000_hpge_metadata" in rules
    assert "plot_geom_rendering" in rules
    assert "plot_geom_hpge_mass" not in rules


def test_optmap_patch_rule_gated_on_config(tmp_path):
    """`optical_maps.lar_patch` inserts `patch_optical_map` ahead of the opt tier.

    Without the key the rule must not exist at all: it would otherwise have no
    input, and productions that do not patch their map must be unaffected.
    """
    steps = ["vtx", "stp", "par", "opt"]

    rules = dag_rule_names(
        default_config, overrides(tmp_path / "nopatch", make_steps=steps)
    )
    assert "build_tier_opt" in rules
    assert "patch_optical_map" not in rules

    # the DAG never reads the maps, so the dummy map stands in for the patch too
    dummy_map = "$_/inputs/simprod/l200cfg01-optmap-dummy.lh5"
    cfg = overrides(tmp_path / "patch", make_steps=steps)
    cfg["paths"]["optical_maps"] = {
        "lar": dummy_map,
        "lar_patch": {"lar_inside": dummy_map},
    }

    rules = dag_rule_names(default_config, cfg)
    assert {"build_tier_opt", "patch_optical_map"} <= rules


# Rerun tests: which jobs Snakemake runs after a change to the metadata. The
# production is run for real, but every job other than the detector-info caches
# only touches its outputs. The working directory is the throwaway one, so the
# job records (in `.snakemake/`) start empty.

PAR_PER_DETECTOR_RULES = {
    "build_hpge_drift_time_map",
    "build_hpge_pulse_shape_library",
    "build_superpulses_from_data",
    "extract_current_pulse_model",
    "extract_electronics_model_pars",
    "convolve_hpge_ideal_pulse_shape_lib",
}
MERGE_RULES = {
    "merge_hpge_drift_time_maps",
    "merge_current_pulse_model_pars",
    "merge_electronics_model_pars",
    "merge_hpge_realistic_psls",
}
CACHE_RULES = {"cache_modelable_hpges", "cache_detector_usabilities"}


def run_prod(
    prod: Mapping, rebuild_caches: bool = False, separately: bool = True
) -> dict[str, list[dict]]:
    """Run *prod* and return the wildcards of the jobs that ran, by rule.

    With *rebuild_caches*, the detector-info caches are rebuilt as after a
    metadata update: in a first Snakemake call if *separately*, otherwise in the
    same call.
    """
    prod["ran"].clear()
    run = partial(
        run_snakemake,
        default_config,
        prod["cfg"],
        workdir=prod["base"],
        executor="local",
    )
    forcerun = CACHE_RULES if rebuild_caches else ()
    if rebuild_caches and separately:
        detinfo = Path(prod["cfg"]["paths"]["pars"]) / "detinfo"
        run(
            targets=[detinfo / "is_modelable.yaml", detinfo / "usability.yaml"],
            forcerun=forcerun,
        )
        forcerun = ()
    run(forcerun=forcerun)
    return dict(prod["ran"])


@pytest.fixture
def recorded_prod(tmp_path, monkeypatch):
    """A production in which every job has run."""
    ran = defaultdict(list)
    run_single_job = local.Executor.run_single_job

    def touch_outputs(job):
        for f in map(Path, [*job.output, *([job.benchmark] if job.benchmark else [])]):
            f.parent.mkdir(parents=True, exist_ok=True)
            f.touch()

    def run_cache_rules_only(executor, job):
        ran[job.rule.name].append(dict(job.wildcards))
        if job.rule.name in CACHE_RULES:
            return run_single_job(executor, job)
        return executor.pool.submit(touch_outputs, job)

    monkeypatch.setattr(local.Executor, "run_single_job", run_cache_rules_only)

    (tmp_path / "workflow").symlink_to(dummyprod / "workflow")
    cfg = overrides(
        tmp_path,
        experiment="legend",
        settings_by_tier={"hit": {"simulate_psd_with_psl": True}},
    )
    prod = {
        "base": tmp_path,
        "cfg": cfg,
        "metadata": Path(cfg["paths"]["metadata"]),
        "ran": ran,
    }
    jobs = run_prod(prod)
    assert {"build_superpulses_from_data", "build_tier_evt"} <= jobs.keys()
    assert run_prod(prod) == {}
    return prod


def edit_status(prod: Mapping, hpge: str, **fields) -> None:
    """Edit the status of *hpge* in the p02 physics runs."""
    status = prod["metadata"] / "datasets/statuses/l200-p02-r006-T%-all-config.yaml"
    data = yaml.safe_load(status.read_text())
    data[hpge] |= fields
    status.write_text(yaml.safe_dump(data))


def test_hpge_cache_loaded_once(recorded_prod, monkeypatch):
    loads = Counter()
    load_dict = dbetto.utils.load_dict

    def counting_load_dict(path, *args, **kwargs):
        loads[Path(path).name] += 1
        return load_dict(path, *args, **kwargs)

    monkeypatch.setattr(dbetto.utils, "load_dict", counting_load_dict)
    run_prod(recorded_prod)
    assert loads["is_modelable.yaml"] == 1
    assert loads["operational_voltage_in_V.yaml"] == 1


def test_usability_change_reruns_tiers(recorded_prod):
    """B99000A is not modelable, before and after."""
    edit_status(recorded_prod, "B99000A", usability="ac")
    # metadata changes are not tracked
    assert run_prod(recorded_prod) == {}

    ran = run_prod(recorded_prod, rebuild_caches=True)
    assert ran.keys() >= CACHE_RULES | {"build_tier_opt", "build_tier_hit"}
    assert "build_tier_evt" in ran
    assert ran.keys().isdisjoint(PAR_PER_DETECTOR_RULES | MERGE_RULES)


def test_same_detector_info_reruns_nothing(recorded_prod):
    """A status field that is not cached changes."""
    edit_status(recorded_prod, "V99000A", processable=False)
    assert run_prod(recorded_prod, rebuild_caches=True).keys() == CACHE_RULES


@pytest.mark.parametrize("separately", [True, False])
def test_modelability_change_reruns_affected_jobs(recorded_prod, config, separately):
    """V99000A is skipped in r001 only.

    When the caches are rebuilt in the same Snakemake call, the jobs downstream
    of them are scheduled before the rebuild and run even if their input does
    not change.
    """
    skip = recorded_prod["metadata"] / "simprod/config/pars/legend/geds/skip"
    (skip / "l200-p02-r001-T%-all-skip.yaml").write_text("V99000A: test\n")
    validity = yaml.safe_load((skip / "validity.yaml").read_text())
    validity[1:1] = [
        {"valid_from": "20220202T000000Z", "apply": ["l200-p02-r001-T%-all-skip.yaml"]},
        {
            "valid_from": "20220302T000000Z",
            "mode": "reset",
            "apply": ["l200-p02-r000-T%-all-skip.yaml"],
        },
    ]
    (skip / "validity.yaml").write_text(yaml.safe_dump(validity))

    ran = run_prod(recorded_prod, rebuild_caches=True, separately=separately)
    assert ran.keys() >= CACHE_RULES
    assert ran["merge_hpge_drift_time_maps"] == [{"runid": "l200-p02-r001-phy"}]
    for rule in MERGE_RULES & ran.keys():
        assert ran[rule] == [{"runid": "l200-p02-r001-phy"}]
    assert "aggregate_hpge_ssd_modeling_info" in ran
    assert ran.keys().isdisjoint(PAR_PER_DETECTOR_RULES)

    affected = {
        simid
        for simid in agg.gen_list_of_all_simids(config)
        if "l200-p02-r001-phy" in agg.get_runlist(config, simid)
    }
    hit_simids = {wc["simid"] for wc in ran["build_tier_hit"]}
    if separately:
        assert "build_tier_opt" not in ran
        assert hit_simids == affected
    else:
        assert hit_simids >= affected


def test_runlist_change_reruns_caches(recorded_prod):
    """A simulation gets a run that is not in the caches yet."""
    simconfig = (
        recorded_prod["metadata"] / "simprod/config/tier/hit/legend/simconfig.yaml"
    )
    data = yaml.safe_load(simconfig.read_text())
    data["phbr_surface_Ra228_to_Ac228"]["runlist"].append("l200-p03-r000-phy")
    simconfig.write_text(yaml.safe_dump(data))

    ran = run_prod(recorded_prod)
    assert ran.keys() >= CACHE_RULES
    assert {"runid": "l200-p03-r000-phy"} in ran["merge_hpge_drift_time_maps"]
