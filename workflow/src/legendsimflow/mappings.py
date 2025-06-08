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

# ruff: noqa: F821, T201
from __future__ import annotations

from legendmeta import LegendMetadata

# TODO: add support for SiPMs


def l200a_mageid_to_detname(mageid: int, on: str | datetime | None = None):
    """Convert MaGe identifier to LEGEND HPGe detector name."""
    lmeta = LegendMetadata()
    chmap = lmeta.channelmap(on=on)

    mageid = int(mageid)
    position = mageid % 100
    string = (mageid % 10000 - position) / 100

    return (
        chmap.map("system", unique=False)
        .geds.map("location.string", unique=False)[string]
        .map("location.position")[position]
        .name
    )


def l200a_detname_to_mageid(detname: str, on: str | datetime | None = None):
    """Convert LEGEND HPGe detector name to MaGe identifier."""
    lmeta = LegendMetadata()
    chmap = lmeta.channelmap(on=on)

    return int(
        1010000
        + 100 * chmap[detname].location.string
        + chmap[detname].location.position
    )
