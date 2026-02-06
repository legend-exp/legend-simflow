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
