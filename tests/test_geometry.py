from __future__ import annotations

import yaml
from dbetto import AttrsDict

from legendsimflow import geometry
from legendsimflow.geometry import DEFAULT_VIS_SCENE


def _cfg(config_dir, experiment="legend"):
    return AttrsDict({"paths": {"config": str(config_dir)}, "experiment": experiment})


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
