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

- event categories: LAr veto, A/E and energy in the FCCD layer. The
  functions return one mask per category, of the events (LAr veto) or of the
  HPGe hits (A/E, FCCD). Each event or hit is in exactly one category of each
  level.
- the energy of the HPGe hits for other dead-layer models.
- the output struct of the pdf tier, and reading the sum of its categories.
- spectra on a grid of parameters, in the layout Refit.jl reads as a
  generalised histogram.
"""

from __future__ import annotations

import functools
import itertools
import logging
import operator
from collections.abc import Callable, Iterable, Mapping
from fnmatch import fnmatch
from pathlib import Path

import awkward as ak
import h5py
import hist
import lh5
import numpy as np
from lgdo import Histogram, Scalar, Struct

from .metadata import get_sanitized_fccd

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
        values[name] = grid_axis_values(first, last, step)
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


def grid_axis_values(first: float, last: float, step: float) -> np.ndarray:
    """Return the values of a grid axis of :func:`build_pdf_grid`."""
    return first + step * np.arange(round((last - first) / step) + 1)


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


def category_levels(
    has_spms_coinc: bool, aoe_model: str | None, fccd_thr: float | None
) -> list[dict[str, str]]:
    """Return the event categories (name -> description), one dict per level.

    The names are the keys of the masks of :func:`lar_categories`,
    :func:`aoe_categories` and :func:`fccd_categories`. A level is left out when
    its data is missing.

    Parameters
    ----------
    has_spms_coinc
        whether the LAr veto flag ``coincident.spms`` is available.
    aoe_model
        ``pulse_lib``, ``single_temp`` or `None` without simulated PSD.
    fccd_thr
        energy threshold in keV of the FCCD tag, or `None` without it.
    """
    levels = []
    if has_spms_coinc:
        levels.append(
            {
                "is_spms_coincident": "light seen by the SiPMs (fails the LAr veto)",
                "not_spms_coincident": "no light seen by the SiPMs (passes the LAr veto)",
            }
        )
    if aoe_model is not None:
        aoe_level = {}
        # the single template has no high-side cut
        if aoe_model == "pulse_lib":
            aoe_level["is_high_aoe"] = "A/E above the high-side cut"
        aoe_level["is_ss_aoe"] = "single-site A/E (passes the A/E cuts)"
        aoe_level["is_low_aoe"] = "A/E below the low-side cut (multi-site)"
        aoe_level["is_invalid_aoe"] = "no valid or simulated A/E"
        levels.append(aoe_level)
    if fccd_thr is not None:
        levels.append(
            {
                "is_fccd_tagged": f"more than {fccd_thr} keV deposited in the FCCD layer",
                "not_fccd_tagged": f"at most {fccd_thr} keV deposited in the FCCD layer",
            }
        )
    return levels


def lar_categories(evt: ak.Array) -> dict[str, np.ndarray]:
    """Return the event masks of the LAr veto categories.

    The categories are defined by ``coincident.spms``, the LAr veto flag of
    the `evt` tier, which is `True` when the SiPMs see enough light:

    - ``is_spms_coincident``: light seen, the event fails the LAr veto;
    - ``not_spms_coincident``: no light seen, the event passes the LAr veto.

    Parameters
    ----------
    evt
        events of the `evt` tier.
    """
    spms = ak.to_numpy(evt.coincident.spms)
    return {"is_spms_coincident": spms, "not_spms_coincident": ~spms}


def aoe_categories(psd: ak.Array, aoe_model: str) -> dict[str, np.ndarray]:
    """Return the hit masks of the A/E categories.

    The flags come from ``geds.psd`` and its ``<aoe_model>`` subtable. The A/E
    classifier is the corrected A/E in units of its resolution, compared in
    the `hit` tier with the low- and high-side A/E cuts of the detector. In
    order of precedence:

    - ``is_invalid_aoe``: no A/E available, i.e. ``psd.is_good`` (valid PSD in
      data) or ``has_aoe`` (A/E simulated) is `False`;
    - ``is_high_aoe``: classifier above the high-side cut (``is_high_aoe``).
      Only with the pulse-shape library, the single template has no
      high-side cut;
    - ``is_ss_aoe``: single-site, classifier above the low-side cut
      (``is_single_site``);
    - ``is_low_aoe``: multi-site, classifier below the low-side cut.

    Each hit is in exactly one category.

    Parameters
    ----------
    psd
        ``geds.psd`` of the `evt` tier. The masks have one entry per hit, in
        the order of ``ak.flatten(geds.energy)``.
    aoe_model
        ``pulse_lib`` (pulse-shape library) or ``single_temp`` (single
        template).
    """

    def per_hit(field):
        return ak.to_numpy(ak.flatten(field, axis=None))

    aoe = psd[aoe_model]
    valid = per_hit(psd.is_good) & per_hit(aoe.has_aoe)
    single_site = per_hit(aoe.is_single_site)

    categories = {}
    high = np.zeros(len(valid), dtype=bool)
    if aoe_model == "pulse_lib":
        high = valid & per_hit(aoe.is_high_aoe)
        categories["is_high_aoe"] = high
    categories["is_ss_aoe"] = valid & ~high & single_site
    categories["is_low_aoe"] = valid & ~high & ~single_site
    categories["is_invalid_aoe"] = ~valid
    return categories


def fccd_categories(
    steps: Mapping[str, np.ndarray], fccd_mm: np.ndarray, energy_thr_kev: float
) -> dict[str, np.ndarray]:
    """Return the hit masks of the FCCD categories.

    The energy in the FCCD layer is the sum of the energy depositions near the
    n+ surface that are closer to it than the FCCD of the detector. These are
    the energies as deposited, before the dead-layer model:

    - ``is_fccd_tagged``: more than `energy_thr_kev` in the FCCD layer;
    - ``not_fccd_tagged``: at most `energy_thr_kev` in the FCCD layer.

    Parameters
    ----------
    steps
        the depositions near the n+ surface, from :func:`flatten_surface_steps`.
    fccd_mm
        FCCD in mm of each detector, in the order of the UIDs given to
        :func:`flatten_surface_steps`.
    energy_thr_kev
        energy threshold in keV.
    """
    in_fccd = steps["distance"] < fccd_mm[steps["det"]]
    edep_fccd = np.bincount(
        steps["hit"],
        weights=steps["edep"] * in_fccd,
        minlength=len(steps["energy_bulk"]),
    )
    tagged = edep_fccd > energy_thr_kev
    return {"is_fccd_tagged": tagged, "not_fccd_tagged": ~tagged}


def fccd_by_uid(
    metadata: Mapping,
    name_to_uid: Mapping[str, int],
    surface_depth: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the sorted detector UIDs and their FCCD in mm.

    The FCCD is read from the detector metadata with
    :func:`.metadata.get_sanitized_fccd`. It is limited to `surface_depth`,
    the depth of the steps stored in ``geds.surface``, with a warning. The
    FCCD of the detector of each hit is then
    ``fccd_mm[np.searchsorted(uids, uid)]``.

    Parameters
    ----------
    metadata
        LEGEND metadata database.
    name_to_uid
        map of the detector names to their UIDs, e.g. ``detector_uids`` of
        the `cvt` tier.
    surface_depth
        ``surface_depth_in_mm`` of the `hit` tier settings.
    """
    names = sorted(name_to_uid, key=name_to_uid.get)
    uids = np.array([name_to_uid[n] for n in names], dtype=np.int64)
    fccd_mm = np.zeros(len(names))
    for i, name in enumerate(names):
        fccd_mm[i] = get_sanitized_fccd(metadata, name)
        # the steps deeper than surface_depth are not stored
        if fccd_mm[i] > surface_depth:
            log.warning(
                "FCCD of %s (%s mm) exceeds surface_depth_in_mm, using %s mm",
                name,
                fccd_mm[i],
                surface_depth,
            )
            fccd_mm[i] = surface_depth
    return uids, fccd_mm


def piecewise_linear_activeness(
    distance: np.ndarray, fccd_mm: np.ndarray | float, dlf: float
) -> np.ndarray:
    """Return the charge-collection efficiency of energy depositions near the n+ surface.

    Same model as :func:`reboost.math.functions.piecewise_linear_activeness`,
    but with one FCCD per deposition: 0 up to ``dlf * fccd_mm`` from the
    surface, rising linearly to 1 at `fccd_mm`, and 1 beyond.

    Parameters
    ----------
    distance
        distance of each deposition to the n+ surface in mm.
    fccd_mm
        FCCD in mm, for each deposition or for all.
    dlf
        dead layer fraction, the fraction of the FCCD that is fully inactive.
    """
    dead = fccd_mm * dlf
    with np.errstate(divide="ignore", invalid="ignore"):
        act = np.clip((distance - dead) / (fccd_mm - dead), 0, 1)
    # 0/0 only at the inner edge of a dead layer as thick as the FCCD
    return np.nan_to_num(act, nan=0.0, copy=False)


def flatten_surface_steps(geds: ak.Array, uids: np.ndarray) -> dict[str, np.ndarray]:
    """Return the HPGe depositions near the n+ surface as flat arrays.

    The hits are numbered in the order of ``ak.flatten(geds.energy)``.

    Parameters
    ----------
    geds
        HPGe hits of the `evt` tier, with the ``surface`` subtable.
    uids
        sorted detector UIDs.

    Returns
    -------
    ``energy_bulk`` (keV), one entry per hit, and one entry per deposition:
    ``hit`` (the number of its hit), ``det`` (the index of its detector in
    `uids`), ``edep`` (keV) and ``distance`` (to the n+ surface, mm).
    """
    surface = geds.surface
    n_steps = ak.to_numpy(ak.flatten(ak.num(surface.edep, axis=-1)))
    hit = np.repeat(np.arange(len(n_steps)), n_steps)
    hit_det = np.searchsorted(uids, ak.to_numpy(ak.flatten(geds.rawid)))
    return {
        "energy_bulk": ak.to_numpy(ak.flatten(surface.energy_bulk)),
        "hit": hit,
        "det": hit_det[hit],
        "edep": ak.to_numpy(ak.flatten(surface.edep, axis=None)),
        "distance": ak.to_numpy(ak.flatten(surface.distance_to_nplus, axis=None)),
    }


def collected_energy(
    steps: Mapping[str, np.ndarray], fccd_mm: np.ndarray, dlf: float
) -> np.ndarray:
    """Return the energy (keV) collected from each HPGe hit, before smearing.

    The collected energy is ``energy_bulk`` plus the depositions near the n+
    surface weighted by their activeness (see
    :func:`piecewise_linear_activeness`). For another dead-layer model, the
    energy of a hit is its nominal (smeared) energy plus the change of its
    collected energy, so the smearing of each hit is kept.

    Parameters
    ----------
    steps
        the depositions near the n+ surface, from :func:`flatten_surface_steps`.
    fccd_mm
        FCCD in mm of each detector, in the order of the UIDs given to
        :func:`flatten_surface_steps`.
    dlf
        dead layer fraction.
    """
    act = piecewise_linear_activeness(steps["distance"], fccd_mm[steps["det"]], dlf)
    return steps["energy_bulk"] + np.bincount(
        steps["hit"], weights=steps["edep"] * act, minlength=len(steps["energy_bulk"])
    )


def format_pdf_struct(
    histograms: Mapping[str, Mapping[str, hist.Hist]],
    descriptions: Mapping[str, str],
) -> Struct:
    """Nest the histograms of the pdf tier by the levels of their path.

    Parameters
    ----------
    histograms
        histograms by path (e.g. ``mul1/not_spms_coincident/is_ss_aoe``) and
        detector group. The groups are the leaves of the output.
    descriptions
        description of each path, stored in the ``description`` attribute of
        its histograms.
    """
    output = Struct()
    for path, per_group in histograms.items():
        for g, h in per_group.items():
            output.add_field(
                f"{path}/{g}", Histogram(h, attrs={"description": descriptions[path]})
            )
    return output


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

    paths = []
    with h5py.File(lh5_file, "r") as f:
        f["pdf"].visititems(
            lambda name, obj: paths.append(name)
            if isinstance(obj, h5py.Group) and "binning" in obj
            else None
        )

    selected = []
    for path in paths:
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
