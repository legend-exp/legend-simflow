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
        # pDz is the full length of a pyg4ometry Tubs
        length = outer_ms.pDz

        # type conversions from pyg4ometry types
        if not isinstance(r_max, float | int):
            r_max = r_max.eval()

        if not isinstance(length, float | int):
            length = length.eval()

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
            f"/RMG/Generator/Confinement/Geometrical/Cylinder/Height {length} mm"
        )

    return lines


def get_minishroud_outer_surface_confine_commands(
    reg: pyg4ometry.geant4.Registry,
    pattern: str | Sequence[str] = "minishroud_tube*",
    surface_max_intersections: int = 10,
) -> list[str]:
    """Extract the commands to sample on the outer surface of the minishrouds (NMS).

    Vertices are sampled on the full surface of the NMS volumes. The ones on the
    inner side (facing the HPGe string) are rejected with one excluded cylinder
    per NMS, whose radius and bottom face lie halfway through the NMS wall. The
    result is the outer side wall and the outer bottom face of each NMS.

    Parameters
    ----------
    reg
        The registry describing the geometry.
    pattern
        The pattern used to search for physical volumes of minishrouds. Each
        solid must be a vertical cylinder closed at the bottom: a ``Tubs``
        minus a narrower ``Tubs`` shifted upwards.
    surface_max_intersections
        Maximum number of intersections of a line with the NMS surface.

    Returns
    -------
    A list of confinement commands for remage.

    """
    string_list = _get_matching_volumes(list(reg.physicalVolumeDict.keys()), pattern)

    if len(string_list) == 0:
        msg = f"no physical volumes matching pattern {pattern} found in the GDML!"
        raise ValueError(msg)

    lines = [
        "/RMG/Generator/Confine Volume",
        "/RMG/Generator/Confinement/SamplingMode SubtractGeometrical",
        "/RMG/Generator/Confinement/FirstSamplingVolume Physical",
        "/RMG/Generator/Confinement/SampleOnSurface true",
        f"/RMG/Generator/Confinement/SurfaceSampleMaxIntersections {surface_max_intersections}",
    ]
    lines += [f"/RMG/Generator/Confinement/Physical/AddVolume {s}" for s in string_list]

    def _val(x):
        return float(x.eval()) if hasattr(x, "eval") else float(x)

    for s in string_list:
        vol = reg.physicalVolumeDict[s]
        solid = vol.logicalVolume.solid

        if not (
            isinstance(solid, pyg4ometry.geant4.solid.Subtraction)
            and isinstance(solid.obj1, pyg4ometry.geant4.solid.Tubs)
            and isinstance(solid.obj2, pyg4ometry.geant4.solid.Tubs)
        ):
            msg = f"solid of physical volume '{s}' is not a Tubs minus a Tubs"
            raise ValueError(msg)

        if any(_val(a) != 0 for a in [*vol.rotation.eval(), *solid.tra2[0].eval()]):
            msg = f"physical volume '{s}' or its inner Tubs is rotated"
            raise ValueError(msg)

        center = [_val(c) for c in vol.position.eval()]
        r_out = _val(solid.obj1.pRMax)
        r_in = _val(solid.obj2.pRMax)
        # pDz is the full length of a pyg4ometry Tubs
        length = _val(solid.obj1.pDz)
        dz_in = _val(solid.tra2[1].eval()[2])

        # z of the outer and inner bottom faces, relative to the NMS center
        z_bottom_out = -length / 2
        z_bottom_in = dz_in - _val(solid.obj2.pDz) / 2
        z_low = (z_bottom_out + z_bottom_in) / 2
        # extend 1 mm above the NMS top
        z_high = length / 2 + 1

        lines += [
            "/RMG/Generator/Confinement/Geometrical/AddExcludedSolid Cylinder",
            f"/RMG/Generator/Confinement/Geometrical/CenterPositionX {center[0]} mm",
            f"/RMG/Generator/Confinement/Geometrical/CenterPositionY {center[1]} mm",
            f"/RMG/Generator/Confinement/Geometrical/CenterPositionZ {center[2] + (z_low + z_high) / 2} mm",
            f"/RMG/Generator/Confinement/Geometrical/Cylinder/OuterRadius {(r_out + r_in) / 2} mm",
            f"/RMG/Generator/Confinement/Geometrical/Cylinder/Height {z_high - z_low} mm",
        ]

    return lines
