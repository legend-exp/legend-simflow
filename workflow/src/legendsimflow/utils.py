# Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
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

from __future__ import annotations

import fnmatch
from datetime import datetime
from pathlib import Path

import legenddataflowscripts as ldfs
import pyg4ometry as pg4
from dbetto import AttrsDict
from legendmeta import LegendMetadata

from . import SimflowConfig


def init_simflow_context(raw_config: dict, workflow) -> AttrsDict:
    if not raw_config:
        msg = "you must set a config file with --configfile"
        raise RuntimeError(msg)

    raw_config.setdefault("benchmark", {"enabled": False})
    ldfs.workflow.utils.subst_vars_in_snakemake_config(workflow, raw_config)
    config = AttrsDict(raw_config)

    # convert all strings in the "paths" block to pathlib.Path
    def _make_path(d):
        for k, v in d.items():
            if isinstance(v, str):
                d[k] = Path(v)
            else:
                d[k] = _make_path(v)
        return d

    config["paths"] = _make_path(config.paths)

    # NOTE: this will attempt a clone of legend-metadata, if the directory does not exist
    # NOTE: don't use lazy=True, we need a fully functional TextDB
    metadata = LegendMetadata(config.paths.metadata)
    if "legend_metadata_version" in config:
        metadata.checkout(config.legend_metadata_version)
    config["metadata"] = metadata

    return AttrsDict(
        {
            "config": config,
            "basedir": workflow.basedir,
            "proctime": datetime.now().strftime("%Y%m%dT%H%M%SZ"),
        }
    )


def _get_matching_volumes(volume_list: list, patterns: str | list) -> list[int]:
    """Get the list of volumes from the GDML. The string can include wildcards."""

    wildcard_list = [patterns] if isinstance(patterns, str) else patterns

    # find all volumes matching at least one pattern
    matched_list = []
    for key in volume_list:
        for name in wildcard_list:
            if fnmatch.fnmatch(key, name):
                matched_list.append(key)
    return matched_list


def get_lar_minishroud_confine_commands(
    reg: pg4.geant4.Registry,
    pattern: str = "minishroud_side*",
    inside: bool = True,
    lar_name: str = "LAr",
) -> list[str]:
    """Extract the commands for the LAr confinement inside/ outside the NMS from the GDML.

    Parameters
    ----------
    reg
        The registry describing the geometry.
    pattern
        The pattern used to search for physical volumes of minishrouds.
    inside
        Where to generate points only inside NMS to exclude them.
    lar_name
        The name of the physical volume of the LAr.

    Returns
    -------
    a list of confinement commands for remage.
    """
    string_list = _get_matching_volumes(list(reg.physicalVolumeDict.keys()), pattern)
    # correct sampling mode
    mode = "IntersectPhysicalWithGeometrical" if inside else "SubtractGeometrical"

    # physical volume sampling
    lines = [f"/RMG/Generator/Confinement/SamplingMode {mode}"]
    lines += [f"/RMG/Generator/Confinement/Physical/AddVolume {lar_name}"]

    for s in string_list:
        vol = reg.physicalVolumeDict[s]

        center = vol.position.eval()
        solid = vol.logicalVolume.solid

        outer_ms = solid.obj1
        r_max = outer_ms.pRMax
        dz = outer_ms.pDz

        # type conversions from pyg4ometry types
        if not isinstance(r_max, float):
            r_max = r_max.eval()

        if not isinstance(dz, float):
            dz = dz.eval()

        command = "AddSolid" if inside else "AddExcludeSolid"
        lines.append(f"/RMG/Generator/Confinement/Geometrical/{command} Cylinder ")

        lines.append(
            f"/RMG/Generator/Confinement/Geometrical/CenterPositionX {center[0]} mm"
        )
        lines.append(
            f"/RMG/Generator/Confinement/Geometrical/CenterPositionY {center[1]} mm"
        )
        lines.append(
            f"/RMG/Generator/Confinement/Geometrical/CenterPositionZ {center[2]} mm"
        )
        lines.append(
            f"/RMG/Generator/Confinement/Geometrical/Cylinder/OuterRadius {r_max} mm"
        )
        lines.append(f"/RMG/Generator/Confinement/Geometrical/Cylinder/Height {dz} mm")

    return lines


def setup_logdir_link(config: SimflowConfig, proctime):
    logdir = Path(config.paths.log)
    logdir.mkdir(parents=True, exist_ok=True)

    # create a handy link to access latest log directory
    link = logdir / "latest"
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(proctime, target_is_directory=True)
