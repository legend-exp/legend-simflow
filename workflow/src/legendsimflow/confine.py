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
from collections.abc import Iterable, Sequence

import pyg4ometry
import pygeomtools


def _get_matching_volumes(
    volume_list: Iterable[str], patterns: str | Sequence[str]
) -> list[str]:
    """Return volumes from `volume_list` whose names match `patterns`.

    Wildcard patterns are supported via :func:`fnmatch.fnmatch`.

    Parameters
    ----------
    volume_list
        List of volume names to search.
    patterns
        Single wildcard pattern string or a list of patterns.

    """
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
    pattern: str | Sequence[str] = "minishroud_tube*",
    inside: bool = True,
    lar_name: str = "liquid_argon",
    outer_radius_in_mm: float | None = None,
    outer_height_in_mm: float | None = None,
) -> list[str]:
    """Extract the commands for the LAr confinement inside/outside the NMS from the GDML.

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
    outer_radius_in_mm
        If provided, gives an outer radius for the confinement. Only supported
        for outside confinement (inside=False).
    outer_height_in_mm
        If provided, gives an outer height for the confinement. Only supported
        for outside confinement (inside=False).

    Returns
    -------
    A list of confinement commands for remage.

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

    if (outer_radius_in_mm is None) != (outer_height_in_mm is None):
        msg = "outer_radius_in_mm and outer_height_in_mm must be provided together"
        raise ValueError(msg)

    if outer_radius_in_mm is not None and outer_height_in_mm is not None:
        if inside:
            msg = (
                "outer_radius and outer_height parameters are only supported for "
                f"outside confinement (inside=False), but inside={inside} was given."
            )
            raise ValueError(msg)
        lines += [
            "/RMG/Generator/Confinement/Geometrical/AddSolid Cylinder",
            f"/RMG/Generator/Confinement/Geometrical/Cylinder/OuterRadius {outer_radius_in_mm} mm",
            f"/RMG/Generator/Confinement/Geometrical/Cylinder/Height {outer_height_in_mm} mm",
        ]

    for s in string_list:
        vol = reg.physicalVolumeDict[s]

        center = vol.position.eval()
        solid = vol.logicalVolume.solid

        # Validate expected geometry structure before accessing attributes
        if not hasattr(solid, "obj1") or solid.obj1 is None:
            msg = (
                f"Expected solid for physical volume '{s}' to have an 'obj1' "
                "attribute representing the outer minishroud cylinder, but it was missing or None."
            )
            raise ValueError(msg)

        outer_ms = solid.obj1

        if not isinstance(outer_ms, pyg4ometry.geant4.solid.Tubs):
            msg = f"Expected solid for physical volume '{s}'.obj1 to be a Tubs,"
            raise ValueError(msg)

        r_max = outer_ms.pRMax
        dz = outer_ms.pDz

        # type conversions from pyg4ometry types
        if not isinstance(r_max, float | int):
            r_max = r_max.eval()

        if not isinstance(dz, float | int):
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


def get_hpge_bulk_confine_commands(reg: pyg4ometry.geant4.Registry) -> list[str]:
    """Remage commands to confine primaries in the bulk of all HPGe detectors.

    Lists every germanium sensitive physical volume registered in the geometry
    (see ``pygeomtools.detectors.get_all_sensvols``) as a remage volume
    confinement. Meant to be used through the ``~function:`` confinement
    mechanism of :func:`legendsimflow.commands.make_remage_macro`.

    Parameters
    ----------
    reg
        The registry describing the geometry.

    Returns
    -------
    list[str]
        Remage confinement commands, one ``AddVolume`` per HPGe physical
        volume.
    """
    volumes = sorted(pygeomtools.detectors.get_all_sensvols(reg, "germanium"))
    if len(volumes) == 0:
        msg = "no germanium sensitive volumes registered in the geometry"
        raise ValueError(msg)

    return ["/RMG/Generator/Confine Volume"] + [
        f"/RMG/Generator/Confinement/Physical/AddVolume {v}" for v in volumes
    ]
