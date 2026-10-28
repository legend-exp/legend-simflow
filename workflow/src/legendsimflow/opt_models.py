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
"""Light models of scintillators, for the ``opt`` tier.

Each model gives the probability that a scintillation photon emitted at a
position produces a photoelectron in any SiPM channel, and replaces the
optical map of a light source (see ``model`` in the ``opt`` tier settings).
The models give the light of all SiPMs together only. Models named
``l1000_*`` hold the layout of the LEGEND-1000 geometry of
legend-pygeom-l1000, and do not apply to other geometries.
"""

from __future__ import annotations

import functools
from typing import TYPE_CHECKING

import numpy as np
import pygeomtools
from numpy.typing import ArrayLike, NDArray

if TYPE_CHECKING:
    from pyg4ometry.geant4 import Registry


def uniform(
    xloc: ArrayLike,
    yloc: ArrayLike,  # noqa: ARG001
    zloc: ArrayLike,  # noqa: ARG001
    *,
    prob: float,
) -> NDArray:
    """Same detection probability `prob` everywhere, e.g. for the PEN plates."""
    return np.full(len(xloc), prob, dtype=np.float64)


@functools.cache
def _l1000_fiber_length(geometry: Registry) -> float:
    """Length (mm) of the straight fiber cores of legend-pygeom-l1000."""
    # the default model builds a fiber module as one tube segment, the
    # detailed model one box per fiber
    lengths = set()
    sensvols = pygeomtools.detectors.get_all_sensvols(geometry, "scintillator")
    for name, det in sensvols.items():
        if det.ntuple_name != "fiber":
            continue
        solid = geometry.physicalVolumeDict[name].logicalVolume.solid
        length_param = {"Tubs": "pDz", "Box": "pZ"}.get(solid.type)
        if length_param is None:
            msg = (
                f"fiber core {name} is a {solid.type}, not a straight fiber "
                "(Tubs or Box) of legend-pygeom-l1000"
            )
            raise ValueError(msg)
        lengths.add(round(solid.evaluateParameterWithUnits(length_param), 6))
    if len(lengths) != 1:
        msg = f"expected fiber cores of a single length, found lengths {sorted(lengths)} mm"
        raise ValueError(msg)
    return lengths.pop()


def l1000_fibers(
    xloc: ArrayLike,  # noqa: ARG001
    yloc: ArrayLike,  # noqa: ARG001
    zloc: ArrayLike,
    *,
    amplitudes: ArrayLike,
    absorption_lengths_in_m: ArrayLike,
    geometry: Registry,
) -> NDArray:
    r"""Detection probability of the light of the LEGEND-1000 fiber cores.

    The light of a fiber is read out by the SiPM channels at its two ends.
    The probability to be detected at one end decreases with the distance
    :math:`d` along the fiber to that end, as a sum of exponentials:

    .. math::

        p(d) = \sum_i a_i \, e^{-d / \lambda_i}

    with :math:`a_i` the `amplitudes` and :math:`\lambda_i` the
    `absorption_lengths_in_m`. The probability of a step is the sum over the
    two ends. The fibers are vertical and straight: they run from the upper
    SiPM holding structure of legend-pygeom-l1000 down, over the length of
    the fiber cores (ntuple ``fiber``) in the `geometry`.
    """
    from pygeoml1000 import hpge_strings  # noqa: PLC0415

    z_top = hpge_strings.z_pos_dict["sipm_upper_holding_structure_upper_end"]
    z_bottom = z_top - _l1000_fiber_length(geometry)
    zloc = np.asarray(zloc)
    return sum(
        a
        * (
            np.exp(-(z_top - zloc) / (1000 * length))
            + np.exp(-(zloc - z_bottom) / (1000 * length))
        )
        for a, length in zip(amplitudes, absorption_lengths_in_m, strict=True)
    )
