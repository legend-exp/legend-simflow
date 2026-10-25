# Copyright (C) 2026 Manuel Huber <info@manuelhu.de>
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
"""Light models of the LEGEND-1000 scintillators, for the ``opt`` tier.

Each model gives the probability that a scintillation photon emitted at a
position produces a photoelectron in a SiPM channel, and replaces the optical
map of a light source (see ``model`` in the ``opt`` tier settings).

The models take the array layout from the configuration of the geometry
(``geom_config``), built with legend-pygeom-l1000: strings, each surrounded by
a fiber barrel of 3 modules with a SiPM channel at the top and at the bottom
end. The channel of module ``m`` of string ``s`` at end ``e`` (``T`` or
``B``) is named ``S<ss><mm><e>``, e.g. ``S0102T``.
"""

from __future__ import annotations

import functools
import re
from types import SimpleNamespace

import numpy as np
from dbetto import AttrsDict
from dbetto.utils import load_dict
from numpy.typing import ArrayLike, NDArray


@functools.cache
def _l1000_strings(geom_config: str) -> tuple[NDArray, NDArray, NDArray]:
    """Number, centre (x, y in mm) and rotation (rad) of the LEGEND-1000 strings."""
    from pygeoml1000 import hpge_strings  # noqa: PLC0415
    from pygeoml1000.config import (  # noqa: PLC0415
        generate_special_metadata,
        load_raw_config,
    )
    from pygeomtools.utils import load_dict_from_config  # noqa: PLC0415

    # as in pygeoml1000.config.resolve_config(), without the channel map, that
    # needs the metadata
    config = load_dict(geom_config)
    meta = SimpleNamespace(
        special_metadata=AttrsDict(
            load_dict_from_config(
                config,
                "special_metadata",
                lambda: generate_special_metadata(
                    load_raw_config(config.get("raw_config"))
                ),
            )
        )
    )
    ids = sorted(meta.special_metadata.hpge_string, key=int)
    strings = [hpge_strings.calculate_string_rotation(i, meta) for i in ids]
    return (
        np.array(ids, dtype=int),
        np.array([[s["x_pos"], s["y_pos"]] for s in strings]),
        np.array([s["string_rot"] for s in strings]),
    )


def _nearest_string(xloc: ArrayLike, yloc: ArrayLike, centres: NDArray) -> NDArray:
    """Index of the string whose axis (`centres`) is closest to each position (mm)."""
    nearest = np.zeros(len(xloc), dtype=int)
    best = np.full(len(xloc), np.inf)
    # a loop instead of a (steps x strings) matrix, to keep the memory low
    for i, (x, y) in enumerate(centres):
        dist = (xloc - x) ** 2 + (yloc - y) ** 2
        nearest[dist < best] = i
        best = np.minimum(best, dist)
    return nearest


def _parse_sipm(sipm: str) -> tuple[int, int, str]:
    """String number, module number and end (``T`` or ``B``) of a SiPM channel."""
    match = re.fullmatch(r"S(\d{2})(\d{2})([TB])", sipm)
    if match is None:
        msg = f"{sipm} is not a LEGEND-1000 SiPM channel name, like S0102T"
        raise ValueError(msg)
    return int(match[1]), int(match[2]), match[3]


def pen_plates(
    xloc: ArrayLike,
    yloc: ArrayLike,
    zloc: ArrayLike,
    *,
    sipm: str,
    prob: float,
    geom_config: str,
) -> NDArray:
    """Detection probability of the light of the PEN plates.

    The light of a plate reaches only the 6 SiPM channels of its string (3
    modules, 2 ends), each with the same probability `prob`. The string of a
    plate is the one with its axis closest to the plate. With `sipm` ``all``,
    the probability is the sum over the 6 channels, ``6 * prob``.
    """
    if sipm == "all":
        return np.full(len(zloc), 6 * prob)

    string, _, _ = _parse_sipm(sipm)
    ids, centres, _ = _l1000_strings(geom_config)
    return np.where(ids[_nearest_string(xloc, yloc, centres)] == string, prob, 0.0)


def fibers(
    xloc: ArrayLike,
    yloc: ArrayLike,
    zloc: ArrayLike,
    *,
    sipm: str,
    amplitudes: ArrayLike,
    absorption_lengths_in_m: ArrayLike,
    fiber_length_in_m: float,
    geom_config: str,
) -> NDArray:
    r"""Detection probability of the light of the fiber cores.

    The light of a fiber reaches only the SiPM channels at the two ends of
    its module. The probability to be detected at one end decreases with the
    distance :math:`d` along the fiber to that end, as a sum of exponentials:

    .. math::

        p(d) = \sum_i a_i \, e^{-d / \lambda_i}

    with :math:`a_i` the `amplitudes` and :math:`\lambda_i` the
    `absorption_lengths_in_m`. The straight fibers run from the upper SiPM
    holding structure of legend-pygeom-l1000 down, over `fiber_length_in_m`
    (1.349 m in legend-pygeom-l1000). The fiber of a step is in the string
    whose axis is closest to it, and in the module whose azimuthal sector
    around that axis holds it. With `sipm` ``all``, the probability is the sum
    over the two ends.
    """
    from pygeoml1000 import hpge_strings  # noqa: PLC0415

    z_top = hpge_strings.z_pos_dict["sipm_upper_holding_structure_upper_end"]
    z_bottom = z_top - 1000 * fiber_length_in_m

    def end_prob(dist_in_mm: NDArray) -> NDArray:
        return sum(
            a * np.exp(-dist_in_mm / (1000 * length))
            for a, length in zip(amplitudes, absorption_lengths_in_m, strict=True)
        )

    if sipm == "all":
        return end_prob(z_top - zloc) + end_prob(zloc - z_bottom)

    string, module, end = _parse_sipm(sipm)
    ids, centres, rots = _l1000_strings(geom_config)
    nearest = _nearest_string(xloc, yloc, centres)
    # module m covers the azimuth range [6 deg, 114 deg] + 120 deg * (m - 1)
    # from the string rotation, the rest is holding structure
    phi = np.arctan2(yloc - centres[nearest, 1], xloc - centres[nearest, 0])
    in_module = np.floor(np.mod(phi - rots[nearest], 2 * np.pi) / (2 * np.pi / 3)) + 1

    prob = end_prob(z_top - zloc) if end == "T" else end_prob(zloc - z_bottom)
    return np.where((ids[nearest] == string) & (in_module == module), prob, 0.0)
