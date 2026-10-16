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
    pdf: Callable[[dict[str, float]], hist.Hist | Histogram],
    axes: Mapping[str, tuple[float, float, float]],
) -> Struct:
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
        value}``, that returns the histogram at this grid point. All histograms
        must have the same binning and hold raw counts, from the same number of
        simulated primaries.
    axes
        for each parameter, in the order of the grid axes, its ``(first, last,
        step)`` values.

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

    grid: dict = {}
    edges = None
    for idx in itertools.product(*(range(len(v)) for v in values.values())):
        h = pdf(
            {name: float(values[name][i]) for name, i in zip(values, idx, strict=True)}
        )
        if not isinstance(h, Histogram):
            h = Histogram(h)

        point_edges = [ax.edges for ax in h.binning]
        if edges is None:
            edges = point_edges
        elif len(point_edges) != len(edges) or not all(
            np.array_equal(a, b) for a, b in zip(point_edges, edges, strict=False)
        ):
            msg = f"the histogram at grid point {idx} has a different binning"
            raise ValueError(msg)

        node = grid
        keys = [f"{name}_{i}" for name, i in zip(values, idx, strict=True)]
        for key in keys[:-1]:
            node = node.setdefault(key, {})
        node[keys[-1]] = h

    def _to_struct(node: dict) -> Struct:
        return Struct(
            {k: _to_struct(v) if isinstance(v, dict) else v for k, v in node.items()}
        )

    grid_info = Struct(
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
    return Struct({"grid_info": grid_info, "grid": _to_struct(grid)})
