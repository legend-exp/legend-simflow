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

import awkward as ak
import numpy as np


def ak_isin(
    elements,
    test_elements,
    *,
    assume_unique=False,
):
    elements_layout = ak.to_layout(elements)

    # match NumPy: treat test_elements as a 1D collection
    test_flat = ak.to_numpy(ak.ravel(test_elements))

    def transformer(layout, *, backend, **kwargs):  # noqa: ARG001
        if layout.is_numpy:
            # layout is a low-level NumpyArray node; convert to NumPy for np.isin
            data = ak.to_numpy(layout)
            out = np.isin(
                data,
                test_flat,
                assume_unique=assume_unique,
            )
            return ak.contents.NumpyArray(
                backend.nplike.asarray(out),
                parameters=layout.parameters,
                backend=backend,
            )
        return None

    return ak.transform(
        transformer,
        elements_layout,
        return_value="simplified",
    )


def select_geds_hit(evt: ak.Array, rawid: int) -> ak.Array:
    """Select the events with a hit in one HPGe detector and keep only that hit.

    The per-hit fields of ``geds`` (for example ``energy`` and the steps in
    ``surface``) lose the hit level: ``geds.energy`` gets one value per event,
    ``geds.surface.edep`` one list of steps per event. The fields with one value
    per event (for example ``geds.energy_sum``) and all the other fields of
    `evt` are kept as they are, for the selected events.

    Parameters
    ----------
    evt
        events of the ``evt`` tier.
    rawid
        the detector, as in the ``geds.rawid`` field.

    Examples
    --------
    >>> hit = select_geds_hit(evt, 1104000)
    >>> hit.geds.energy  # one energy per event, in keV
    >>> hit.geds.surface.distance_to_nplus  # one list of steps per event, in mm
    >>> hit.coincident.spms  # one value per selected event
    """
    is_det = evt.geds.rawid == rawid
    if ak.any(ak.sum(is_det, axis=-1) > 1):
        msg = f"detector {rawid} has more than one hit in an event"
        raise ValueError(msg)

    in_det = ak.any(is_det, axis=-1)
    evt = evt[in_det]
    is_det = is_det[in_det]

    def _select(field):
        if field.fields:
            return ak.zip({k: _select(field[k]) for k in field.fields}, depth_limit=1)
        if field.ndim == 1:
            return field
        return ak.flatten(field[is_det], axis=1)

    # packed: repeated computations on the steps (e.g. a PDF grid scan) run faster
    return ak.with_field(evt, ak.to_packed(_select(evt.geds)), "geds")
