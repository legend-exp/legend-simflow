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
"""Spectra on a grid of parameters, in the LH5 format of a generalised histogram."""

from __future__ import annotations

import itertools
from collections.abc import Callable, Mapping

import hist
import numpy as np
from lgdo import Histogram, Scalar, Struct


def build_pdf_grid(
    pdf: Callable[
        [dict[str, float]], hist.Hist | Histogram | Mapping[str, hist.Hist | Histogram]
    ],
    axes: Mapping[str, tuple[float, float, float]],
) -> Struct | dict[str, Struct]:
    """Evaluate a PDF at every point of an evenly spaced grid of parameters.

    The result is a generalised histogram as read by Refit.jl: a ``grid_info``
    struct with the name and the ``first``, ``last`` and ``step`` values of each
    parameter (``axis_0``, ``axis_1``, ...), and a ``grid`` struct with the
    histogram of each grid point at ``grid/<name0>_<i>/<name1>_<j>/...``,
    nested in the order of `axes`, with indices counted from 0. Write it with
    ``lh5.write``.

    Parameters
    ----------
    pdf
        function of the parameter values, given as a dictionary ``{name:
        value}``, that returns the histogram at this grid point, or a dictionary
        of histograms ``{key: histogram}`` to build several PDFs in one scan,
        with the same keys at every grid point. The histograms of each PDF must
        have the same binning at all grid points and hold raw counts, from the
        same number of simulated primaries.
    axes
        for each parameter, in the order of the grid axes, its ``(first, last,
        step)`` values.

    Returns
    -------
    the generalised histogram, or a dictionary ``{key: generalised histogram}``
    if `pdf` returns a dictionary.

    Examples
    --------
    >>> grid = build_pdf_grid(pdf, {"fccd_sigma": (-5, 5, 1), "dlf": (0, 1, 0.25)})
    >>> lh5.write(grid, "pdf_grid/mul1/all", "grid.lh5")
    """
    values = {}
    for name, (first, last, step) in axes.items():
        n = round((last - first) / step) + 1
        values[name] = first + step * np.arange(n)
        if n < 2 or not np.isclose(values[name][-1], last):
            msg = f"axis {name} does not go from {first} to {last} in steps of {step}"
            raise ValueError(msg)

    # the grid tree and the bin edges of each PDF
    grids: dict = {}
    edges: dict = {}
    for idx in itertools.product(*(range(len(v)) for v in values.values())):
        result = pdf(
            {name: float(values[name][i]) for name, i in zip(values, idx, strict=True)}
        )
        is_single = not isinstance(result, Mapping)
        hists = {None: result} if is_single else result

        if grids and set(hists) != set(grids):
            msg = f"the PDFs at grid point {idx} differ from those at the first point"
            raise ValueError(msg)

        point = [f"{name}_{i}" for name, i in zip(values, idx, strict=True)]
        for key, point_hist in hists.items():
            h = (
                point_hist
                if isinstance(point_hist, Histogram)
                else Histogram(point_hist)
            )

            point_edges = [ax.edges for ax in h.binning]
            if key not in edges:
                edges[key] = point_edges
            elif len(point_edges) != len(edges[key]) or not all(
                np.array_equal(a, b)
                for a, b in zip(point_edges, edges[key], strict=False)
            ):
                msg = f"the histogram {key} at grid point {idx} has a different binning"
                raise ValueError(msg)

            node = grids.setdefault(key, {})
            for level in point[:-1]:
                node = node.setdefault(level, {})
            node[point[-1]] = h

    def _to_struct(node: dict) -> Struct:
        return Struct(
            {k: _to_struct(v) if isinstance(v, dict) else v for k, v in node.items()}
        )

    def _grid_info() -> Struct:
        return Struct(
            {
                f"axis_{d}": Struct(
                    {
                        "name": Scalar(name),
                        "values": Struct(
                            {
                                k: Scalar(float(v))
                                for k, v in zip(
                                    ("first", "last", "step"), axes[name], strict=True
                                )
                            }
                        ),
                    }
                )
                for d, name in enumerate(values)
            }
        )

    out = {
        key: Struct({"grid_info": _grid_info(), "grid": _to_struct(tree)})
        for key, tree in grids.items()
    }
    return out[None] if is_single else out
