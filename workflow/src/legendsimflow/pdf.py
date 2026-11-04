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
"""Probability density functions (PDFs) of the pdf tier.

- the output struct of the pdf tier, and reading the sum of its categories.
  The categories are defined in :mod:`legendsimflow.selections`.
- spectra on a grid of parameters, in the layout Refit.jl reads as a
  generalised histogram.
"""

from __future__ import annotations

import functools
import itertools
import logging
import operator
import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from fnmatch import fnmatch
from pathlib import Path

import h5py
import hist
import lh5
import numpy as np
from lgdo import Histogram, Scalar, Struct

log = logging.getLogger(__name__)


def build_pdf_grid(
    pdf: Callable[
        [dict[str, float]], hist.Hist | Histogram | Mapping[str, hist.Hist | Histogram]
    ],
    axes: Mapping[str, tuple[float, float, float]],
) -> Struct | dict[str, Struct]:
    """Evaluate a PDF at every point of an evenly spaced grid of parameters.

    The result is an LGDO struct that Refit.jl reads as a generalised histogram:
    a ``grid_info`` struct with the name and the ``first``, ``last`` and ``step`` values of each
    parameter (``axis_0``, ``axis_1``, ...), and a ``grid`` struct with the
    histogram of each grid point at ``grid/<name0>_<i>/<name1>_<j>/...``,
    nested in the order of `axes`, with indices counted from 0. Write it with
    ``lh5.write``.

    Parameters
    ----------
    pdf
        function of the parameter values, given as a dictionary ``{name:
        value}``, that returns the histogram at this grid point, or a dictionary
        of histograms ``{key: histogram}`` to build several PDFs in one scan.
        The histograms of a PDF must have the same binning at all grid points
        and hold raw counts, from the same number of simulated primaries.
    axes
        for each parameter, in the order of the grid axes, its ``(first, last,
        step)`` values.

    Returns
    -------
    the grid, or a dictionary ``{key: grid}`` if `pdf` returns a dictionary.

    Examples
    --------
    >>> grid = build_pdf_grid(pdf, {"fccd_sigma": (-5, 5, 1), "dlf": (0, 1, 0.25)})
    >>> lh5.write(grid, "pdf_grid/mul1/all", "grid.lh5")
    """
    grid_info = Struct()
    values = {}
    for d, (name, (first, last, step)) in enumerate(axes.items()):
        values[name] = first + step * np.arange(round((last - first) / step) + 1)
        grid_info.add_field(f"axis_{d}/name", Scalar(name))
        for key, value in (("first", first), ("last", last), ("step", step)):
            grid_info.add_field(f"axis_{d}/values/{key}", Scalar(float(value)))

    grids = {}
    for idx in itertools.product(*(range(len(v)) for v in values.values())):
        hists = pdf(
            {name: float(values[name][i]) for name, i in zip(values, idx, strict=True)}
        )
        if not isinstance(hists, Mapping):
            hists = {None: hists}

        point = "/".join(f"{name}_{i}" for name, i in zip(values, idx, strict=True))
        for key, h in hists.items():
            if key not in grids:
                grids[key] = Struct({"grid_info": grid_info, "grid": Struct()})
            grids[key]["grid"].add_field(
                point, h if isinstance(h, Histogram) else Histogram(h)
            )

    for grid in grids.values():
        validate_grid(grid)

    # one grid if pdf returns one histogram
    return grids.get(None, grids)


def validate_grid(grid: Struct) -> None:
    """Check a grid as returned by :func:`build_pdf_grid`.

    Raises a :class:`ValueError` if an axis does not go from ``first`` to
    ``last`` in whole steps, if a grid point has no histogram, or if the
    histograms have different binnings.

    Parameters
    ----------
    grid
        the grid, built or read from a file with ``lh5.read``.
    """
    n_values = {}
    for axis in grid["grid_info"].values():
        name = axis["name"].value
        name = name.decode() if isinstance(name, bytes) else name  # bytes when read
        first, last, step = (axis["values"][k].value for k in ("first", "last", "step"))
        n_values[name] = round((last - first) / step) + 1
        if n_values[name] < 2 or not np.isclose(
            first + (n_values[name] - 1) * step, last
        ):
            msg = f"axis {name} does not go from {first} to {last} in steps of {step}"
            raise ValueError(msg)

    edges = None
    for idx in itertools.product(*(range(n) for n in n_values.values())):
        point = "/".join(f"{name}_{i}" for name, i in zip(n_values, idx, strict=True))
        if point not in grid["grid"]:
            msg = f"grid point {point} has no histogram"
            raise ValueError(msg)

        point_edges = [ax.edges for ax in grid["grid"][point].binning]
        if edges is None:
            edges = point_edges
        elif len(point_edges) != len(edges) or not all(
            np.array_equal(a, b) for a, b in zip(point_edges, edges, strict=False)
        ):
            msg = f"the histogram at grid point {point} has a different binning"
            raise ValueError(msg)


def group_detectors(
    groups: Mapping[str, str] | None, detectors: Mapping[str, int]
) -> dict[str, np.ndarray]:
    """Return the IDs of the detectors in each group, plus the implicit ``all`` group.

    Parameters
    ----------
    groups
        the ``detector_groups`` setting of the pdf tier: group name to regular
        expression, matched against the detector names with
        :func:`re.fullmatch`.
    detectors
        detector name to ID: the UID in simulations, the rawid in data.
    """
    groups = dict(groups or {})
    if groups.get("all", ".*") != ".*":
        log.warning(
            "detector_groups sets 'all' to %r; it will be overridden by the implicit all-detector group '.*'",
            groups["all"],
        )
    groups["all"] = ".*"

    group_ids = {}
    for name, pattern in groups.items():
        rgx = re.compile(pattern)
        matched = sorted(i for n, i in detectors.items() if rgx.fullmatch(n))
        if not matched:
            log.warning("detector group %r matches no detectors", name)
        group_ids[name] = np.asarray(matched, dtype=np.int64)
    return group_ids


def new_histograms(
    levels: Sequence[Mapping[str, str]], groups: Iterable[str]
) -> tuple[dict[tuple[str, ...], dict[str, hist.Hist]], hist.Hist]:
    """Return the empty histograms of the pdf tier.

    One energy histogram, 1 keV bins from 0 to 6000 keV, per combination of
    categories (one of each level of `levels`) and detector group, keyed by
    the combination and the group. And the 2D histogram of the
    multiplicity-2 events, with the same binning on both axes.

    Parameters
    ----------
    levels
        the category levels, see :func:`.selections.category_levels`.
    groups
        the detector group names.
    """
    # boost-histogram only provides Double (float64) storage
    histograms = {
        combo: {g: hist.new.Reg(6000, 0, 6000).Double() for g in groups}
        for combo in itertools.product(*levels)
    }
    mul2_hist = hist.new.Reg(6000, 0, 6000).Reg(6000, 0, 6000).Double()
    return histograms, mul2_hist


def fill_histograms(
    histograms: Mapping[tuple[str, ...], Mapping[str, hist.Hist]],
    energy: np.ndarray,
    detector: np.ndarray,
    in_category: Mapping[str, np.ndarray],
    group_ids: Mapping[str, np.ndarray],
) -> None:
    """Fill the histograms of :func:`new_histograms` with multiplicity-1 events.

    Parameters
    ----------
    histograms
        the histograms, by combination of categories and group.
    energy
        energy in keV of the hit of each event.
    detector
        ID of the detector of the hit of each event.
    in_category
        event mask of each category, see :mod:`legendsimflow.selections`.
    group_ids
        IDs of the detectors of each group, see :func:`group_detectors`.
    """
    in_group = {g: np.isin(detector, ids) for g, ids in group_ids.items()}
    for combo, per_group in histograms.items():
        sel = np.ones(len(energy), dtype=bool)
        for c in combo:
            sel &= in_category[c]
        for g, h in per_group.items():
            h.fill(energy[sel & in_group[g]])


def format_pdf_struct(
    histograms: Mapping[tuple[str, ...], Mapping[str, hist.Hist]],
    levels: Sequence[Mapping[str, str]],
    mul2_hist: hist.Hist | None = None,
) -> Struct:
    """Nest the histograms of the pdf tier by category and detector group.

    The histogram of a combination of categories and a group goes to
    ``mul1/<category>/.../<group>``, with a ``description`` attribute listing
    the categories. `mul2_hist` goes to ``mul2/all``.

    Parameters
    ----------
    histograms
        the histograms of :func:`new_histograms`, by combination of
        categories and group.
    levels
        the category levels, name to description, see
        :func:`.selections.category_levels`.
    mul2_hist
        the 2D histogram of the multiplicity-2 events.
    """
    descriptions = {k: v for level in levels for k, v in level.items()}
    output = Struct()
    for combo, per_group in histograms.items():
        desc = "multiplicity 1" + "".join(f", {descriptions[c]}" for c in combo)
        for g, h in per_group.items():
            output.add_field(
                "/".join(("mul1", *combo, g)),
                Histogram(h, attrs={"description": desc}),
            )
    if mul2_hist is not None:
        output.add_field(
            "mul2/all",
            Histogram(
                mul2_hist,
                attrs={
                    "description": "multiplicity-2 events: (E_low, E_high) for pairs of ON detectors"
                },
            ),
        )
    return output


def histogram_paths(lh5_file: str | Path) -> list[str]:
    """Return the paths of the histograms below ``pdf/`` in a pdf tier file."""
    paths = []
    with h5py.File(lh5_file, "r") as f:
        f["pdf"].visititems(
            lambda name, obj: paths.append(name)
            if isinstance(obj, h5py.Group) and "binning" in obj
            else None
        )
    return paths


def merge_pdf_files(files: Iterable[str | Path], output: str | Path) -> None:
    """Sum the histograms of several pdf tier files into one file.

    The histograms at the same path are added, so the files must hold the
    same paths and binning: the PDFs of one simulation over several
    productions, or the data histograms of several runs. The root scalars
    ``nr_sim_events`` and ``livetime_in_s`` are summed too, when present.

    Parameters
    ----------
    files
        the pdf tier files to sum.
    output
        the file to write.
    """
    files = [str(f) for f in files]
    pdf = Struct()
    for path in histogram_paths(files[0]):
        h = lh5.read(f"pdf/{path}", files[0])
        for f in files[1:]:
            h.weights.nda += lh5.read(f"pdf/{path}", f).weights.nda
        pdf.add_field(path, h)
    lh5.write(pdf, "pdf", output, wo_mode="w")

    for name in ("nr_sim_events", "livetime_in_s"):
        if name in lh5.ls(files[0]):
            total = sum(lh5.read(name, f).value for f in files)
            lh5.write(Scalar(total), name, output, wo_mode="append")


def read_pdf(
    selection: str | Iterable[str],
    lh5_file: str | Path,
    group: str = "all",
) -> hist.Hist:
    """Read the sum of the 1D histograms of the pdf tier below one or more paths.

    The categories of the pdf tier do not overlap, so any selection is a sum
    of categories. A path in `selection` selects every category below it,
    and ``*`` matches any category at its level.

    Parameters
    ----------
    selection
        path or paths below ``pdf/``, without the detector group.
    lh5_file
        output file of the pdf tier.
    group
        detector group.

    Examples
    --------
    >>> read_pdf("mul1", "pdf.lh5")  # all multiplicity-1 events
    >>> read_pdf("mul1/not_spms_coincident", "pdf.lh5")  # after the LAr veto
    >>> read_pdf("mul1/*/is_ss_aoe", "pdf.lh5")  # after the A/E cuts only
    >>> read_pdf(
    ...     ["mul1/*/is_ss_aoe", "mul1/*/is_low_aoe"],  # after the high-side A/E cut
    ...     "pdf.lh5",
    ...     group="icpc",
    ... )
    """
    patterns = [selection] if isinstance(selection, str) else list(selection)
    patterns = [pat.strip("/").split("/") for pat in patterns]

    selected = []
    for path in histogram_paths(lh5_file):
        *levels, g = path.split("/")
        if g == group and any(
            len(levels) >= len(pat)
            and all(fnmatch(level, p) for level, p in zip(levels, pat, strict=False))
            for pat in patterns
        ):
            selected.append(path)

    if not selected:
        msg = f"no histogram of group {group!r} matches {selection!r} in {lh5_file}"
        raise ValueError(msg)

    return functools.reduce(
        operator.add,
        (lh5.read_as(f"pdf/{path}", str(lh5_file), "hist") for path in selected),
    )
