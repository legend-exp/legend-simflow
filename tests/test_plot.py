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

import re
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import awkward as ak
import hist
import lh5
import matplotlib.pyplot as plt
import numpy as np
from lgdo import Array, Table
from matplotlib.backends.backend_pdf import PdfPages

from legendsimflow import plot


def test_decorate():
    fig, _ = plt.subplots()
    plot.decorate(fig)
    plt.close(fig)


def test_set_empty():
    fig, ax = plt.subplots()
    plot.set_empty(ax)
    texts = [t.get_text() for t in ax.texts]
    assert "empty!" in texts
    plt.close(fig)


def test_save_page(tmp_path):
    pdf_path = tmp_path / "out.pdf"
    with PdfPages(pdf_path) as pdf:
        plot.save_page(pdf, lambda: plt.subplots()[0])
    assert pdf_path.exists()


def test_plot_hist_nonempty():
    h = hist.Hist(hist.axis.Regular(10, 0, 10))
    h.fill([1, 2, 3, 4, 5])
    fig, ax = plt.subplots()
    result = plot.plot_hist(h, ax)
    assert result is True
    plt.close(fig)


def test_plot_hist_empty():
    h = hist.Hist(hist.axis.Regular(10, 0, 10))
    fig, ax = plt.subplots()
    result = plot.plot_hist(h, ax)
    assert result is False
    texts = [t.get_text() for t in ax.texts]
    assert "empty!" in texts
    plt.close(fig)


def test_plot_hist_with_nans():
    h = hist.Hist(hist.axis.Regular(10, 0, 10))
    h.fill([1, 2, 3])
    fig, ax = plt.subplots()
    result = plot.plot_hist(h, ax, n_nans=2)
    assert result is True
    plt.close(fig)


def test_n_nans_none():
    arr = ak.Array([1.0, 2.0, 3.0])
    assert plot.n_nans(arr) == 0


def test_n_nans_some():
    arr = ak.Array([1.0, float("nan"), 3.0, float("nan")])
    assert plot.n_nans(arr) == 2


def test_read_concat_wempty_no_files():
    result = plot.read_concat_wempty([], "some_table")
    assert result is None


def test_read_concat_wempty_with_table(tmp_path):
    f = tmp_path / "test.lh5"
    table = Table(size=3)
    table.add_field("energy", Array(np.array([1.0, 2.0, 3.0])))
    lh5.write(table, "hit/detectors", f)

    result = plot.read_concat_wempty([f], "detectors")
    assert result is not None
    assert len(result) == 3


def test_read_concat_wempty_missing_table(tmp_path):
    f = tmp_path / "test.lh5"
    table = Table(size=2)
    lh5.write(table, "hit/other_table", f)

    result = plot.read_concat_wempty([f], "detectors")
    assert result is None


def test_all_savefig_calls_are_preceded_by_decorate():
    repo_root = Path(__file__).resolve().parents[1]
    files = [
        repo_root / "workflow/src/legendsimflow/plot.py",
        *sorted((repo_root / "workflow/src/legendsimflow/scripts").rglob("*.py")),
    ]

    decorate_pattern = re.compile(r"\b(?:plot\.)?decorate\(")
    lookback = 8

    for file in files:
        lines = file.read_text().splitlines()

        for lineno, line in enumerate(lines, start=1):
            code_only = line.split("#", 1)[0]
            if "savefig(" not in code_only:
                continue

            prior_lines = [
                ln.split("#", 1)[0]
                for ln in lines[max(0, lineno - lookback - 1) : lineno]
            ]
            assert decorate_pattern.search("\n".join(prior_lines)), (
                f"{file}:{lineno} has savefig call without decorate in previous "
                f"{lookback} lines"
            )
