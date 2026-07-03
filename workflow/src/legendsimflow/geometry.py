# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Helpers producing the `stp` geometry validation plots."""

from __future__ import annotations

import copy
import logging
import os
from collections.abc import Mapping
from pathlib import Path

import dbetto

from . import SimflowConfig, patterns

# NOTE: the heavy plotting/geometry packages (matplotlib, VTK via pyg4ometry /
# pygeomtools, pygeoml200) are imported lazily inside the two functions below,
# not at module level, so the module stays light for `load_vis_scene` and tests.

# default scene passed to :func:`pygeomtools.viewer.visualize`. It hides the large
# opaque enclosures and the outer fiber barrel to expose the detector array,
# makes the inner fiber barrel translucent, and uses a fixed 3/4 view of the
# LEGEND-200 array (world frame, mm). Productions can override it (wholesale, per top-level key) with a
# `<experiment>-vis-config.yaml` file next to the geometry config in the
# metadata; see :func:`load_vis_scene`. The `fine_mesh` key is handled by the
# rendering rule (it must be applied before the geometry is built) and the
# `export_and_exit` output path is injected there too.
DEFAULT_VIS_SCENE: dict = {
    "window_size": [1600, 2600],
    "fine_mesh": True,
    "light": {"pos": [-5000, 0, 5000], "shadow": True},
    "default": {
        "focus": [2.67, -48.24, 711.62],
        "up": [0.26, -0.27, 1.0],
        "camera": [-2440.94, 2563.88, 2184.31],
        "parallel": False,
    },
    "color_overrides": {
        "cryostat_steel": False,
        "liquid_argon": False,
        "gaseous_argon": False,
    },
}


def load_vis_scene(config: SimflowConfig) -> dict:
    """Return the geometry rendering scene for the current experiment.

    Starts from ``DEFAULT_VIS_SCENE`` and merges (shallow, per top-level key) an
    optional per-experiment override read from
    `<paths.config>/geom/<experiment>-vis-config.yaml` in the metadata.
    """
    scene = copy.deepcopy(DEFAULT_VIS_SCENE)
    override = patterns.geom_vis_config_filename(config)
    if override.exists():
        scene |= dbetto.utils.load_dict(override)
    return scene


def make_hpge_mass_plot(
    config: SimflowConfig, geom_config: Mapping, output: str
) -> None:
    """Write the simulated-vs-measured HPGe mass comparison plot to *output*.

    Compares the mass of each detector built by {mod}`pygeoml200` to the measured
    mass in `legend-metadata` (or the public testdata masses for a public
    geometry).
    """
    import matplotlib.pyplot as plt  # noqa: PLC0415
    from pygeoml200.hpge_mass_check import plot_hpge_mass_comparison  # noqa: PLC0415

    from .plot import decorate  # noqa: PLC0415

    # pygeoml200.hpge_mass_check pulls in VTK (via pygeomtools), which switches
    # matplotlib to a Qt backend; force the head-less Agg backend back
    plt.switch_backend("agg")

    os.environ["LEGEND_METADATA"] = str(config.paths.metadata)

    # silence pygeoml200's expected noise (per-detector dummy-enrichment
    # warnings, public-geometry notice) that would otherwise spam the log
    logging.getLogger("pygeoml200").setLevel(logging.ERROR)

    fig, ax = plt.subplots(figsize=(13, 3.5))
    plot_hpge_mass_comparison(geom_config, ax=ax)
    decorate(fig, rotate=True)
    fig.savefig(output, bbox_inches="tight")


def render_geometry(config: SimflowConfig, geom_config: Mapping, output: str) -> None:
    """Render the geometry off-screen to *output* (PNG) using :func:`load_vis_scene`.

    Rebuilds the geometry with the light-weight segmented fiber model (the
    simulation GDML uses the detailed one, far too heavy to render). Rendering
    goes through the software OSMesa backend, so no GPU or X server is needed.
    """
    from pyg4ometry import config as meshconfig  # noqa: PLC0415
    from pygeoml200 import cli, core  # noqa: PLC0415
    from pygeomtools import viewer  # noqa: PLC0415

    # software OSMesa off-screen rendering; must be set before the render window
    # is created
    os.environ["VTK_DEFAULT_OPENGL_WINDOW"] = "vtkOSOpenGLRenderWindow"
    os.environ["LEGEND_METADATA"] = str(config.paths.metadata)

    # silence pygeoml200's expected noise (per-detector dummy-enrichment
    # warnings, public-geometry notice) that would otherwise spam the log
    logging.getLogger("pygeoml200").setLevel(logging.ERROR)

    scene = load_vis_scene(config)
    if scene.pop("fine_mesh", False):  # must be applied before building the geometry
        meshconfig.setGlobalMeshSliceAndStack(100)

    registry = core.construct(
        assemblies=cli._parse_assemblies(geom_config.get("assemblies")),
        use_detailed_fiber_model=False,
        config=geom_config,
        public_geometry=geom_config.get("public_geom", False),
    )

    # `viewer._export_png` refuses to overwrite, so clear a stale target first
    Path(output).unlink(missing_ok=True)
    scene["export_and_exit"] = str(output)
    viewer.visualize(registry, scene)
