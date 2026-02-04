# ruff: noqa: I002

# Copyright (C) 2026 Luigi Pertoldi <gipert@pm.me>,
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

from pathlib import Path

import awkward as ak
import hist
import matplotlib.pyplot as plt
from lgdo import lh5
from matplotlib.backends.backend_pdf import PdfPages

from legendsimflow import nersc, plot

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

hit_files = args.input
output_pdf = args.output[0]
simid = args.wildcards.simid


def fig(table):
    fig = plt.figure(figsize=(12, 6))

    data = plot.read_concat_wempty(hit_files, table)

    if len(data) == 0:
        ax = fig.add_subplot()
        plot.set_empty(ax)
        return fig

    outer = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1, 1])
    gs_top = outer[0].subgridspec(1, 2, width_ratios=[1, 1])
    gs_bot = outer[1].subgridspec(1, 2, width_ratios=[1, 1])

    # time
    ax = fig.add_subplot(gs_top[0, 0])
    h_time = hist.new.Reg(300, 0, 3000, name="photoelectron $t - t_0$ (ns)").Double()
    h_time.fill_flattened(data.time - data.t0)
    plot.plot_hist(h_time, ax)
    ax.set_ylabel("counts / 10 ns")
    ax.set_yscale("log")

    # number of photoelectrons
    ax = fig.add_subplot(gs_top[0, 1])
    h_npe = hist.new.Reg(500, 0, 500, name="number of photoelectrons").Double()
    h_npe.fill(ak.count(data.time, axis=-1))
    plot.plot_hist(h_npe, ax)
    ax.set_ylabel("counts / 1 pe")
    ax.set_yscale("log")

    ax = fig.add_subplot(gs_bot[0, 0])
    plot.set_empty(ax)

    # usability
    ax = fig.add_subplot(gs_bot[0, 1])
    plot.set_empty(ax)
    # vals, counts = np.unique(data.usability, return_counts=True)
    # plt.pie(
    #     counts, labels=[mutils.decode_usability(v) for v in vals], autopct="%1.1f%%"
    # )
    # ax.set_aspect("equal")  # keep it circular

    fig.suptitle(f"{simid}: {table} hits")
    return fig


# prepare a pdf file with a plot per page
tables: list[str] = [
    Path(t).name
    for t in sorted(lh5.ls(hit_files[0], "hit/"))
    if not Path(t).name.startswith("__")
]
fig_builders = [lambda t=t: fig(t) for t in tables]

with PdfPages(output_pdf) as pdf:
    for make_fig in fig_builders:
        plot.save_page(pdf, make_fig)
