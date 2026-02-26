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

import awkward as ak
import hist
import matplotlib.pyplot as plt
from lgdo.lh5 import LH5Iterator

from legendsimflow import nersc, plot

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

BUFFER_LEN = "100*MB"


def _gimme_ehist():
    return hist.new.Reg(500, 0, 5000, name="energy (keV)").Double()


cvt_file = args.input
output_pdf = args.output[0]
simid = args.wildcards.simid


def _fill_ehist(h, evt_chunk, mask):
    evt = evt_chunk.view_as("ak")

    energy = evt.geds.energy

    # apply mask at energy level (mask is jagged: per-hit)
    energy = energy[mask]

    # remove empty arrays
    energy = energy[ak.count(energy, axis=-1) > 0]

    # compute event total energy
    energy = ak.sum(energy, axis=-1)

    h.fill(energy)


def _plot_ehist(ax, mask, **kwargs):
    kwargs = {"fill": True} | kwargs

    h = _gimme_ehist()
    it = LH5Iterator(
        cvt_file,
        "evt",
        buffer_len=BUFFER_LEN,
    )
    for evt_chunk in it:
        _fill_ehist(h, evt_chunk, mask(evt_chunk.view_as("ak")))

    plot.plot_hist(h, ax=ax, **kwargs)


fig = plt.figure(figsize=(14, 12))

outer = fig.add_gridspec(nrows=3, ncols=1, height_ratios=[1, 1, 1])
gs_top = outer[0].subgridspec(1, 1)
gs_mid = outer[1].subgridspec(1, 1)
gs_bot = outer[2].subgridspec(1, 2, width_ratios=[1, 1])

ax = fig.add_subplot(gs_top[0, 0])

_plot_ehist(
    ax,
    lambda evt: evt.coincident.geds,
    color="black",
    fill=False,
    linewidth=1,
    label="evt.coincident.geds",
)
base_mask = lambda evt: evt.coincident.geds & evt.geds.is_good_channel  # noqa: E731
_plot_ehist(
    ax,
    base_mask,
    color="tab:gray",
    fill=False,
    label="geds.is_good_channel",
)
_plot_ehist(
    ax,
    lambda evt: base_mask(evt) & (evt.geds.multiplicity == 1),
    color="silver",
    label="... geds.multiplicity == 1",
)
_plot_ehist(
    ax,
    lambda evt: base_mask(evt) & (evt.geds.multiplicity == 1) & ~evt.coincident.spms,
    color="tab:blue",
    label="... ~coincident.spms",
)

ax.set_yscale("log")
ax.legend()

ax = fig.add_subplot(gs_mid[0, 0])

base_mask = (  # noqa: E731
    lambda evt: evt.coincident.geds
    & (evt.geds.multiplicity == 1)
    & evt.geds.is_good_channel
    & evt.geds.has_aoe
)
_plot_ehist(
    ax,
    base_mask,
    color="silver",
    label="geds.is_good_channel & geds.has_aoe & geds.multiplicity == 1",
)
_plot_ehist(
    ax,
    lambda evt: base_mask(evt) & evt.geds.is_single_site,
    color="tab:green",
    label="... geds.is_single_site",
)
_plot_ehist(
    ax,
    lambda evt: base_mask(evt) & evt.geds.is_single_site & ~evt.coincident.spms,
    color="tab:red",
    label="... ~coincident.spms",
)

ax.set_yscale("log")
ax.legend()

ax = fig.add_subplot(gs_bot[0, 0])
h = hist.new.IntCategory(range(10), name="geds multiplicity").Double()
it = LH5Iterator(
    cvt_file, "evt", buffer_len=BUFFER_LEN, field_mask=["geds/multiplicity"]
)
for evt_chunk in it:
    h.fill(evt_chunk.view_as("ak").geds.multiplicity)
plot.plot_hist(h, ax)

ax = fig.add_subplot(gs_bot[0, 1])
h = hist.new.IntCategory(range(60), name="spms multiplicity").Double()
it = LH5Iterator(
    cvt_file, "evt", buffer_len=BUFFER_LEN, field_mask=["spms/multiplicity"]
)
for evt_chunk in it:
    h.fill(evt_chunk.view_as("ak").spms.multiplicity)
plot.plot_hist(h, ax)

fig.suptitle(f"{simid}: evt tier")

fig.savefig(output_pdf)
