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
from __future__ import annotations

import fnmatch

import pyg4ometry


def _get_matching_volumes(volume_list: list, patterns: str | list) -> list[str]:
    """Get the list of volumes from the GDML. The string can include wildcards."""

    wildcard_list = [patterns] if isinstance(patterns, str) else patterns

    # find all volumes matching at least one pattern
    matched_list = []
    matched_set = set()
    for key in volume_list:
        for name in wildcard_list:
            if fnmatch.fnmatch(key, name):
                if key not in matched_set:
                    matched_list.append(key)
                    matched_set.add(key)
                break
    return matched_list


def get_lar_minishroud_confine_commands(
    reg: pyg4ometry.geant4.Registry,
    pattern: str = "minishroud_tube*",
    inside: bool = True,
    lar_name: str = "liquid_argon",
) -> list[str]:
    """Extract the commands for the LAr confinement inside/ outside the NMS from the GDML.

    Parameters
    ----------
    reg
        The registry describing the geometry.
    pattern
        The pattern used to search for physical volumes of minishrouds.
    inside
        If True, generate points inside the minishroud (NMS) volumes; if False,
        exclude the minishroud volumes from the generation region.
    lar_name
        The name of the physical volume of the LAr.

    Returns
    -------
    a list of confinement commands for remage.
    """
    string_list = _get_matching_volumes(list(reg.physicalVolumeDict.keys()), pattern)

    if len(string_list) == 0:
        msg = f"no physical volumes matching pattern {pattern} found in the GDML!"
        raise ValueError(msg)

    # correct sampling mode
    mode = "IntersectPhysicalWithGeometrical" if inside else "SubtractGeometrical"

    # physical volume sampling

    lines = [
        "/RMG/Generator/Confine Volume",
        f"/RMG/Generator/Confinement/SamplingMode {mode}",
    ]
    lines += [f"/RMG/Generator/Confinement/Physical/AddVolume {lar_name}"]

    for s in string_list:
        vol = reg.physicalVolumeDict[s]

        center = vol.position.eval()
        solid = vol.logicalVolume.solid

        # Validate expected geometry structure before accessing attributes
        if not hasattr(solid, "obj1") or solid.obj1 is None:
            msg = (
                f"Expected solid for physical volume '{s}' to have an 'obj1' "
                "attribute representing the outer minishroud cylinder, but it was  missing or None."
            )
            raise ValueError(msg)

        outer_ms = solid.obj1

        if not isinstance(outer_ms, pyg4ometry.geant4.solid.Tubs):
            msg = f"Expected solid for physical volume '{s}'.obj1 to be a Tubs,"
            raise ValueError(msg)

        r_max = outer_ms.pRMax
        dz = outer_ms.pDz

        # type conversions from pyg4ometry types
        if not isinstance(r_max, float):
            r_max = r_max.eval()

        if not isinstance(dz, float):
            dz = dz.eval()

        command = "AddSolid" if inside else "AddExcludedSolid"
        lines.append(f"/RMG/Generator/Confinement/Geometrical/{command} Cylinder")

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
        lines.append(
            f"/RMG/Generator/Confinement/Geometrical/Cylinder/Height {2 * dz} mm"
        )

    return lines
