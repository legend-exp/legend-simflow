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

from collections import Counter
from pathlib import Path

import awkward as ak
import hist
import lh5
import matplotlib.pyplot as plt
import numpy as np
from lh5 import LH5Iterator
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

    bw = 5
    h_ene = hist.new.Reg(int(5000 / bw), 0, 5000, name="energy (keV)").Double()
    # 0.5 keV bins, the zoom inset is cut from it
    h_ene_fine = hist.new.Reg(10000, 0, 5000, name="energy (keV)").Double()
    h_aoe = hist.new.Reg(200, 0, 2, name="A/E").Double()
    h_dt = hist.new.Reg(
        300, 0, 3000, name="drift time · $t_{max(A)} - t_0$ (ns)"
    ).Double()
    n_rows = n_aoe = n_dt = n_nans_aoe = n_nans_dt = 0
    usability = Counter()
    # hits above 1 MeV for the A/E classifier panels, in file order
    high_e = []
    n_high_e = 0

    # some hit files do not contain the table because there were no hits
    files = [f for f in hit_files if f"hit/{table}" in lh5.ls(f, "hit/")]
    chunks = (
        LH5Iterator(
            files,
            f"hit/{table}",
            field_mask=[
                "energy",
                "usability",
                "psd/single_temp/aoe",
                "psd/single_temp/aoe_corr",
                "psd/single_temp/drift_time_amax",
                "psd/single_temp/is_single_site",
            ],
            buffer_len="100*MB",
        )
        if files
        else []
    )
    for chunk in chunks:
        data = chunk.view_as("ak")
        n_rows += len(data)

        h_ene.fill(data.energy)
        h_ene_fine.fill(data.energy)

        aoe_raw = data.psd.single_temp.aoe_corr[data.energy > 100]
        h_aoe.fill(aoe_raw)
        n_aoe += len(aoe_raw)
        n_nans_aoe += n_nans(aoe_raw)

        dt = data.psd.single_temp.drift_time_amax[data.energy > 100]
        h_dt.fill(dt)
        n_dt += len(dt)
        n_nans_dt += n_nans(dt)

        vals, counts = np.unique(data.usability, return_counts=True)
        usability.update(dict(zip(vals.tolist(), counts.tolist(), strict=True)))

        if n_high_e < 10_000:
            # copy: the iterator overwrites its buffer with the next chunk
            sel = ak.copy(data[data.energy > 1000][: 10_000 - n_high_e])
            high_e.append(sel)
            n_high_e += len(sel)

    if n_rows == 0:
        ax = fig.add_subplot()
        plot.set_empty(ax)
        return fig

    outer = fig.add_gridspec(nrows=3, ncols=1, height_ratios=[1, 1, 1])
    gs_top = outer[0].subgridspec(1, 2, width_ratios=[4, 1])
    gs_mid = outer[1].subgridspec(1, 3, width_ratios=[1, 1, 0.6])
    gs_bot = outer[2].subgridspec(1, 2, width_ratios=[1, 4])

    # energy
    ax = fig.add_subplot(gs_top[0, 0])
    plot.plot_hist(h_ene, ax)
    ax.set_ylabel(f"counts / {bw} keV")
    ax.set_yscale("log")

    # energy inset
    ax = fig.add_subplot(gs_top[0, 1])

    # find tallest gamma peak above 1 MeV
    _h = h_ene[1000j:]
    if _h.sum() != 0:
        argmax_e = _h.axes[0].centers[_h.view().argmax()]

        # +-20 keV around the centre of the 5 keV bin, which is the edge
        # 10 * k + 5 of the 0.5 keV bins
        i = 10 * h_ene.axes[0].index(argmax_e) + 5
        h_ene_z = h_ene_fine[i - 40 : i + 40]
        plot.plot_hist(h_ene_z, ax, flow="none")
    else:
        # no events above 1 MeV; mark inset as empty
        plot.set_empty(ax)

    ax.set_ylabel("counts / 0.5 keV")
    ax.set_yscale("log")

    # A/E
    ax = fig.add_subplot(gs_mid[0, 0])
    if n_aoe > 0:
        plotted = plot.plot_hist(
            h_aoe,
            ax,
            color="tab:red",
            label="energy > 100 keV",
            n_nans=n_nans_aoe,
        )
        if plotted:
            ax.legend()
    else:
        plot.set_empty(ax)

    # drift time
    ax = fig.add_subplot(gs_mid[0, 1])
    if n_dt > 0:
        plotted = plot.plot_hist(
            h_dt,
            ax,
            color="tab:orange",
            label="energy > 100 keV",
            n_nans=n_nans_dt,
        )
        if plotted:
            ax.legend()
    else:
        plot.set_empty(ax)

    # usability
    ax = fig.add_subplot(gs_mid[0, 2])
    vals = sorted(usability)
    labels = [mutils.decode_usability(v) for v in vals]
    plt.pie(
        [usability[v] for v in vals],
        labels=labels,
        colors=[plot.USABILITY_COLOR[lab] for lab in labels],
        autopct="%1.1f%%",
    )
    ax.set_aspect("equal")  # keep it circular

    # A/E classifier
    _d = ak.concatenate(high_e)
    _ss = _d.psd.single_temp.is_single_site

    ax1 = fig.add_subplot(gs_bot[0, 0])
    h_aoec = hist.new.Reg(100, -20, 5).Double().fill(_d.psd.single_temp.aoe[_ss])
    plot.plot_hist(
        hist.new.Reg(100, -20, 5).Double().fill(_d.psd.single_temp.aoe[~_ss]),
        ax1,
        color="tab:gray",
        flow="none",
        orientation="horizontal",
    )
    plotted = plot.plot_hist(
        h_aoec,
        ax1,
        color="tab:red",
        label="energy > 1 MeV",
        flow="none",
        orientation="horizontal",
    )

    def norm(y):
        scale = h_aoec.sum() * h_aoec.axes[0].widths[0]
        return scale * np.exp(-0.5 * y**2) / np.sqrt(2 * np.pi)

    if plotted:
        y = np.linspace(-5, 5, 1000)
        ax1.plot(norm(y), y, label=r"$\mathcal{N}(0,1)$")
        ax1.legend(loc="lower left")

    ax1.invert_xaxis()
    ax1.set_ylabel("A/E classifier")

    ax2 = fig.add_subplot(gs_bot[0, 1], sharey=ax1)
    ax2.scatter(
        _d.energy[_ss],
        _d.psd.single_temp.aoe[_ss],
        s=2,
        color="tab:red",
        label="is_single_site",
    )
    ax2.scatter(_d.energy[~_ss], _d.psd.single_temp.aoe[~_ss], s=2, color="tab:gray")
    ax2.set_xlim(1000, None)
    ax2.set_ylim(-20, 5)
    ax2.set_xlabel("energy (keV)")
    if len(_ss) > 0:
        ax2.legend()
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
