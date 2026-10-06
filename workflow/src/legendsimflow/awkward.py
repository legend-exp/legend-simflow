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

from collections.abc import Callable, Sequence

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


def ak_by_channel(
    values: ak.Array,
    uids: ak.Array,
    channels: Sequence[int],
    reduce: Callable[..., ak.Array] | None = None,
) -> ak.Array:
    """Arrange the hits of each event in a fixed list of channels.

    Parameters
    ----------
    values
        a value per hit, or a list of values per hit, with shape ``events x hits``
        (``events x hits x values``), as read through a time-coincidence map.
    uids
        the channel of each hit, with shape ``events x hits``.
    channels
        the channels of the output, sorted in ascending order. Hits of other
        channels are dropped.
    reduce
        how to combine the scalar values of several hits in the same channel,
        e.g. :func:`awkward.sum`; called with ``axis=1``. If ``None``, `values`
        must hold lists, and the lists of the hits in a channel are joined in
        hit order.

    Returns
    -------
    an array of shape ``events x channels`` (``events x channels x values``).
    Channels without hits hold an empty list, or what `reduce` gives for no
    hits.
    """
    channels = np.asarray(channels)
    n_events, n_channels = len(uids), len(channels)

    # position of each hit's channel in `channels`, and whether it is there
    flat_uids = ak.to_numpy(ak.flatten(uids))
    slot = np.searchsorted(channels, flat_uids)
    known = slot < n_channels
    known[known] = channels[slot[known]] == flat_uids[known]

    # one key per (event, channel); the stable sort keeps the hit order
    event = np.repeat(np.arange(n_events), ak.to_numpy(ak.num(uids)))
    key = (event * n_channels + slot)[known]
    order = np.argsort(key, kind="stable")

    # group the hits by key, with an empty group for every key without hits
    flat = ak.flatten(values, axis=1)[known][order]
    per_slot = ak.unflatten(flat, np.bincount(key, minlength=n_events * n_channels))
    per_slot = (
        ak.flatten(per_slot, axis=2) if reduce is None else reduce(per_slot, axis=1)
    )

    return ak.unflatten(per_slot, np.full(n_events, n_channels))
