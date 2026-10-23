from __future__ import annotations

import copy
from pathlib import Path

import pytest
from dbetto import AttrsDict

from legendsimflow import aggregate, defaults, metadata

docs = Path(__file__).parents[1] / "docs/source/manual"


def test_merge():
    default = {"a": 1, "b": {"c": 2, "d": [1, 2]}, "e": "x"}
    snapshot = copy.deepcopy(default)

    out = defaults.merge({"a": 10, "b": {"d": [3]}, "f": None}, default)
    assert isinstance(out, AttrsDict)
    assert out == {"a": 10, "b": {"c": 2, "d": [3]}, "e": "x", "f": None}

    # a nested block replaces a plain value and vice versa
    assert defaults.merge({"e": {"g": 1}}, default).e == {"g": 1}
    assert defaults.merge({"b": None}, default).b is None

    assert defaults.merge(None, default) == default
    assert defaults.merge({}, {}) == {}

    # the defaults are never modified
    out.b.d.append(4)
    out.b.c = 99
    assert default == snapshot


def test_make_steps_default():
    assert defaults.CONFIG["make_steps"] == aggregate.STEPS_ORDERED


def test_elecmod_scan_settings():
    for k, v in defaults.ELECMOD_SCAN_SETTINGS.items():
        if k not in ("tau_start", "tau_limits", "max_calls"):
            assert v == defaults.PAR_SETTINGS["elecmod"][k]


def _keys(block):
    for k, v in block.items():
        yield k
        if isinstance(v, dict):
            yield from _keys(v)


@pytest.mark.parametrize(
    "block",
    [
        defaults.CONFIG,
        *defaults.TIER_SETTINGS.values(),
        *defaults.PAR_SETTINGS.values(),
    ],
)
def test_defaults_are_documented(block):
    text = "".join(f.read_text() for f in docs.glob("*.md"))
    for key in _keys(block):
        assert f"`{key}`" in text or f"{key}:" in text, key


def test_tier_settings_defaults(config):
    user = config.metadata.simprod.config.tier.hit[config.experiment].settings
    settings = metadata.get_tier_settings(config, "hit")

    for key, value in defaults.TIER_SETTINGS["hit"].items():
        assert settings[key] == user.get(key, value)
    # required fields come from the user settings only
    assert settings.dead_layer_fraction == user.dead_layer_fraction

    assert metadata.get_tier_settings(config, "no-such-tier") == {}


def test_par_settings_defaults(config):
    # no modeling settings file in the test metadata
    settings = metadata.get_par_settings(config, "modeling")
    assert settings == defaults.PAR_SETTINGS["modeling"]

    settings = metadata.get_par_settings(
        config, "elecmod", default=defaults.ELECMOD_SCAN_SETTINGS
    )
    assert settings.max_calls == defaults.ELECMOD_SCAN_SETTINGS["max_calls"]
