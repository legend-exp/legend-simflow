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

from collections.abc import Iterable
from datetime import date
from pathlib import Path

import awkward as ak
import hist
import matplotlib.pyplot as plt
import numpy as np
from lgdo import lh5

USABILITY_COLOR = {
    "off": "tab:red",
    "on": "tab:green",
    "ac": "tab:blue",
}


def decorate(fig):
    fig.text(
        1,
        0,
        f"legend-simflow · {date.today().isoformat()}",
        ha="right",
        va="bottom",
        fontsize=8,
        color="0.4",
    )


def save_page(pdf, make_fig):
    fig = make_fig()
    decorate(fig)
    pdf.savefig(fig)
    plt.close(fig)


def set_empty(ax):
    ax.text(
        0.5,
        0.5,
        "empty!",
        transform=ax.transAxes,
        ha="center",
        va="center",
        color="0.6",
        fontsize=20,
    )


def plot_hist(h, ax, n_nans: int | None = None, **kwargs):
    kwargs = {"flow": "show", "yerr": False} | kwargs

    if h.sum() != 0:
        if n_nans is not None and n_nans > 0:
            h_nans = h.copy(deep=True)
            h_nans.reset()
            h_nans[hist.overflow] = n_nans

            opts = {
                "facecolor": "none",  # optional: transparent fill
                "edgecolor": "tab:red",
                "hatch": "///",
                "flow": "show",
                "yerr": False,
                "label": "np.nan",
            }
            h_nans.plot(ax=ax, **(kwargs | opts))

        h.plot(ax=ax, **kwargs)
    else:
        set_empty(ax)


def read_concat_wempty(files: Iterable[str | Path], table: str) -> ak.Array | None:
    # some hit files might not contain the table because there were no hits
    arrays = []
    for f in files:
        if f"hit/{table}" in lh5.ls(f, "hit/"):
            arrays.append(lh5.read_as(f"hit/{table}", f, library="ak"))

    return ak.concatenate(arrays, axis=0) if arrays != [] else None


def n_nans(array):
    return int(ak.sum(np.isnan(array)))
