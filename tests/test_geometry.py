from __future__ import annotations

import yaml
from dbetto import AttrsDict

from legendsimflow import geometry
from legendsimflow.geometry import DEFAULT_GEOM_EXECUTABLE, DEFAULT_VIS_SCENE


def _cfg(config_dir, experiment="legend"):
    return AttrsDict({"paths": {"config": str(config_dir)}, "experiment": experiment})


def _write_geom_config(config_dir, contents, experiment="legend"):
    geom = config_dir / "geom"
    geom.mkdir(exist_ok=True)
    (geom / f"{experiment}-geom-config.yaml").write_text(yaml.safe_dump(contents))


def test_geom_executable_default(tmp_path):
    """A template config without an `executable` field selects the L200 generator."""
    _write_geom_config(tmp_path, {"public_geom": True})
    assert geometry.geom_executable(_cfg(tmp_path)) == DEFAULT_GEOM_EXECUTABLE


def test_geom_executable_override(tmp_path):
    """The `executable` field selects the experiment's geometry generator."""
    _write_geom_config(
        tmp_path,
        {"executable": "legend-pygeom-l1000"},
        experiment="l1000dsg01",
    )
    cfg = _cfg(tmp_path, experiment="l1000dsg01")
    assert geometry.geom_executable(cfg) == "legend-pygeom-l1000"


def test_load_vis_scene_default(tmp_path):
    """Without a metadata override, the built-in default scene is returned."""
    scene = geometry.load_vis_scene(_cfg(tmp_path))
    assert scene == DEFAULT_VIS_SCENE
    # a deep copy is returned: mutating it must not affect the default
    scene["window_size"][0] = -1
    assert DEFAULT_VIS_SCENE["window_size"][0] != -1


def test_load_vis_scene_override(tmp_path):
    """A per-experiment metadata file overrides the default per top-level key."""
    geom = tmp_path / "geom"
    geom.mkdir()
    (geom / "legend-vis-config.yaml").write_text(
        yaml.safe_dump({"window_size": [10, 20]})
    )

    scene = geometry.load_vis_scene(_cfg(tmp_path))
    assert scene["window_size"] == [10, 20]
    # keys not present in the override are kept from the default
    assert scene["color_overrides"] == DEFAULT_VIS_SCENE["color_overrides"]
