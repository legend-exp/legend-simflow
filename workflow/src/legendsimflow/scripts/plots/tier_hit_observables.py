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

import hist
import matplotlib.pyplot as plt
import numpy as np
from lgdo import lh5
from matplotlib.backends.backend_pdf import PdfPages

from legendsimflow import metadata as mutils
from legendsimflow import nersc, plot
from legendsimflow.plot import n_nans

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

hit_files = args.input
output_pdf = args.output[0]
simid = args.wildcards.simid


def fig(table):
    fig = plt.figure(figsize=(14, 12))

    data = plot.read_concat_wempty(hit_files, table)

    if data is None or len(data) == 0:
        ax = fig.add_subplot()
        plot.set_empty(ax)
        return fig

    outer = fig.add_gridspec(nrows=3, ncols=1, height_ratios=[1, 1, 1])
    gs_top = outer[0].subgridspec(1, 2, width_ratios=[4, 1])
    gs_mid = outer[1].subgridspec(1, 3, width_ratios=[1, 1, 0.6])
    gs_bot = outer[2].subgridspec(1, 2, width_ratios=[1, 4])

    # energy
    ax = fig.add_subplot(gs_top[0, 0])
    bw = 5
    h_ene = hist.new.Reg(int(5000 / bw), 0, 5000, name="energy (keV)").Double()
    h_ene.fill(data.energy)
    plot.plot_hist(h_ene, ax)
    ax.set_ylabel(f"counts / {bw} keV")
    ax.set_yscale("log")

    # energy inset
    ax = fig.add_subplot(gs_top[0, 1])

    # find tallest gamma peak above 1 MeV
    _h = h_ene[1000j:]
    if _h.sum() != 0:
        argmax_e = _h.axes[0].centers[_h.view().argmax()]

        h_ene_z = hist.new.Reg(
            80, argmax_e - 20, argmax_e + 20, name="energy (keV)"
        ).Double()
        h_ene_z.fill(data.energy)
        plot.plot_hist(h_ene_z, ax, flow="none")
    else:
        # no events above 1 MeV; mark inset as empty
        plot.set_empty(ax)

    ax.set_ylabel("Counts / 0.5 keV")
    ax.set_yscale("log")

    # A/E
    ax = fig.add_subplot(gs_mid[0, 0])
    aoe_raw = data.aoe_raw[data.energy > 100]
    h_aoe = hist.new.Reg(200, 0, 2, name="A/E").Double()
    h_aoe.fill(aoe_raw)

    if len(aoe_raw) > 0:
        plotted = plot.plot_hist(
            h_aoe,
            ax,
            color="tab:red",
            label="energy > 100 keV",
            n_nans=n_nans(data.aoe_raw),
        )
        if plotted:
            ax.legend()
    else:
        plot.set_empty(ax)

    # drift time
    ax = fig.add_subplot(gs_mid[0, 1])
    h_dt = hist.new.Reg(
        300, 0, 3000, name="drift time · $t_{max(A)} - t_0$ (ns)"
    ).Double()
    dt = data.drift_time_amax[data.energy > 100]
    h_dt.fill(dt)

    if len(dt) > 0:
        plotted = plot.plot_hist(
            h_dt,
            ax,
            color="tab:orange",
            label="energy > 100 keV",
            n_nans=n_nans(data.drift_time_amax),
        )
        if plotted:
            ax.legend()
    else:
        plot.set_empty(ax)

    # usability
    ax = fig.add_subplot(gs_mid[0, 2])
    vals, counts = np.unique(data.usability, return_counts=True)
    labels = [mutils.decode_usability(v) for v in vals]
    plt.pie(
        counts,
        labels=labels,
        colors=[plot.USABILITY_COLOR[lab] for lab in labels],
        autopct="%1.1f%%",
    )
    ax.set_aspect("equal")  # keep it circular

    # A/E classifier
    _d = data[data.energy > 1000]

    ax1 = fig.add_subplot(gs_bot[0, 0])
    h_aoec = hist.new.Reg(100, -20, 5).Double()
    h_aoec.fill(_d.aoe)
    plot.plot_hist(
        h_aoec,
        ax1,
        color="tab:red",
        label="energy > 1 MeV",
        flow="none",
        orientation="horizontal",
    )

    def norm(y):
        scale = 0.5 * h_aoec.sum() * h_aoec.axes[0].widths[0]
        return scale * np.exp(-0.5 * y**2) / np.sqrt(2 * np.pi)

    y = np.linspace(-5, 5, 1000)
    ax1.plot(norm(y), y, label=r"$\mathcal{N}(0,1)$")

    ax1.invert_xaxis()
    ax1.set_ylabel("A/E classifier")

    ax1.legend()

    ax2 = fig.add_subplot(gs_bot[0, 1], sharey=ax1)
    ax2.scatter(_d.energy, _d.aoe, s=2, color="tab:red")
    ax2.set_xlim(1000, None)
    ax2.set_ylim(-20, 5)
    ax2.set_xlabel("energy (keV)")
    ax2.grid()

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
