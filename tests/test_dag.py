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
from collections.abc import Mapping
from pathlib import Path

import pytest
import yaml
from snakemake import api as smkapi
from snakemake.exceptions import MissingInputException, WorkflowError

dummyprod = Path(__file__).parent / "dummyprod"
default_config = dummyprod / "simflow-config.yaml"
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
    "extract_current_pulse_model",
    "merge_current_pulse_model_pars",
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

    ``settings_by_tier`` maps a tier name to the settings keys to overwrite for
    *experiment* (e.g. ``{"evt": {"skip_opt": True}}``). Editing the copy keeps
    the committed metadata untouched.
    """
    meta = base / "inputs"
    shutil.copytree(dummyprod / "inputs", meta)
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


def dag_rule_names(configfile: Path, config_overrides: Mapping) -> set[str]:
    """Return the set of rule names making up the resolved workflow DAG.

    Built with the touch executor (see the module docstring), so rules both
    upstream and downstream of the modelable-HPGe checkpoint appear, independent
    of what exists on disk. Raises if the DAG cannot be resolved (e.g. an input
    has no producing rule).
    """
    output = smkapi.OutputSettings(verbose=False)
    with smkapi.SnakemakeApi(output) as api:
        wf_api = api.workflow(
            snakefile=dummyprod / "workflow/Snakefile",
            workdir=dummyprod,
            config_settings=smkapi.ConfigSettings(
                configfiles=(configfile,), config=dict(config_overrides)
            ),
            storage_settings=smkapi.StorageSettings(),
            resource_settings=smkapi.ResourceSettings(cores=all_cores),
        )
        dag = wf_api.dag()
        # the touch executor dumps the full job list to stdout; swallow it
        with contextlib.redirect_stdout(io.StringIO()):
            dag.execute_workflow(executor="touch")
        # the public API exposes no DAG accessor, so reach into the workflow
        # object for the resolved graph
        return {job.rule.name for job in wf_api._workflow.dag.jobs}


def test_dag(tmp_path):
    """The default (full) DAG resolves without error and covers every tier."""
    rules = dag_rule_names(default_config, overrides(tmp_path))
    assert "all" in rules
    assert set(TIER_BUILD_RULE.values()) <= rules


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


# the PSD switches gate per-detector rules downstream of the modelable-HPGe
# checkpoint; the touch executor (see dag_rule_names) is what expands them
STEPS_TO_HIT = ["vtx", "stp", "par", "opt", "hit"]


def test_simulate_psd_with_psl_toggles_psl_rules(tmp_path):
    """The hit-tier `simulate_psd_with_psl` setting gates the PSL build chain."""
    on = dag_rule_names(
        l1000_config,
        overrides(
            tmp_path / "on",
            make_steps=STEPS_TO_HIT,
            experiment="l1000dsg01",
            settings_by_tier={"hit": {"simulate_psd_with_psl": True}},
        ),
    )
    off = dag_rule_names(
        l1000_config,
        overrides(
            tmp_path / "off",
            make_steps=STEPS_TO_HIT,
            experiment="l1000dsg01",
            settings_by_tier={"hit": {"simulate_psd_with_psl": False}},
        ),
    )
    assert on >= PSL_PSD_RULES
    assert PSL_PSD_RULES.isdisjoint(off)
    # flipping the switch changes nothing but the PSL rules
    assert on - off == PSL_PSD_RULES


def test_simulate_psd_toggles_dtmap_rules(tmp_path):
    """The hit-tier `simulate_psd` setting gates the drift-time map build chain."""
    on = dag_rule_names(
        l1000_config,
        overrides(
            tmp_path / "on",
            make_steps=STEPS_TO_HIT,
            experiment="l1000dsg01",
            settings_by_tier={"hit": {"simulate_psd": True}},
        ),
    )
    off = dag_rule_names(
        l1000_config,
        overrides(
            tmp_path / "off",
            make_steps=STEPS_TO_HIT,
            experiment="l1000dsg01",
            settings_by_tier={"hit": {"simulate_psd": False}},
        ),
    )
    assert on >= STD_PSD_RULES
    assert STD_PSD_RULES.isdisjoint(off)
    # flipping the switch changes nothing but the drift-time map rules
    assert on - off == STD_PSD_RULES


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
        l1000_config,
        overrides(
            tmp_path / "skip",
            make_steps=steps,
            experiment="l1000dsg01",
            settings_by_tier={"evt": {"skip_opt": True}},
        ),
    )
    assert {"build_tier_evt", "build_tier_hit"} <= rules
    assert "build_tier_opt" not in rules

    with pytest.raises(MissingInputException):
        dag_rule_names(
            l1000_config,
            overrides(
                tmp_path / "noskip",
                make_steps=steps,
                experiment="l1000dsg01",
                settings_by_tier={"evt": {"skip_opt": False}},
            ),
        )


def test_skip_hit_drops_hit_tier(tmp_path):
    """The evt-tier `skip_hit` setting drops the hit tier from the DAG."""
    steps = ["vtx", "stp", "par", "opt", "evt"]  # note: no hit
    rules = dag_rule_names(
        l1000_config,
        overrides(
            tmp_path / "skip",
            make_steps=steps,
            experiment="l1000dsg01",
            settings_by_tier={"evt": {"skip_hit": True}},
        ),
    )
    assert {"build_tier_evt", "build_tier_opt"} <= rules
    assert "build_tier_hit" not in rules

    with pytest.raises(MissingInputException):
        dag_rule_names(
            l1000_config,
            overrides(
                tmp_path / "noskip",
                make_steps=steps,
                experiment="l1000dsg01",
                settings_by_tier={"evt": {"skip_hit": False}},
            ),
        )


def test_skip_opt_and_hit_are_mutually_exclusive(tmp_path):
    """Skipping both the opt and hit tiers is rejected at DAG-build time."""
    with pytest.raises(WorkflowError, match="skip_opt and skip_hit"):
        dag_rule_names(
            l1000_config,
            overrides(
                tmp_path,
                make_steps=["vtx", "stp", "par", "evt"],
                experiment="l1000dsg01",
                settings_by_tier={"evt": {"skip_opt": True, "skip_hit": True}},
            ),
        )


def test_geom_plots_scheduled(tmp_path):
    """The geometry validation plots are scheduled with the stp tier."""
    rules = dag_rule_names(default_config, overrides(tmp_path))
    assert "plot_geom" in rules
