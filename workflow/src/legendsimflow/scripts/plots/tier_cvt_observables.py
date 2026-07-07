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
import lh5
import matplotlib.pyplot as plt
from lh5 import LH5Iterator

from legendsimflow import nersc, plot

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

BUFFER_LEN = "100*MB"


def _gimme_ehist():
    return hist.new.Reg(500, 0, 5000, name="energy (keV)").Double()


cvt_file = args.input
output_pdf = args.output[0]
simid = args.wildcards.simid

has_rc = "evt/spms/rc_energy" in lh5.ls(str(cvt_file), "evt/spms/")


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


fig = plt.figure(figsize=(14, 16))

outer = fig.add_gridspec(nrows=4, ncols=1, height_ratios=[1, 1, 1, 1])
gs_top = outer[0].subgridspec(1, 1)
gs_mid = outer[1].subgridspec(1, 1)
gs_bot = outer[2].subgridspec(1, 2, width_ratios=[1, 1])
gs_pe = outer[3].subgridspec(1, 2, width_ratios=[1, 1])

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
    & evt.geds.psd.single_temp.has_aoe
    & evt.geds.psd.is_good
)
_plot_ehist(
    ax,
    base_mask,
    color="silver",
    label="geds.is_good_channel & geds.psd.single_temp.has_aoe & geds.psd.is_good & geds.multiplicity == 1",
)
_plot_ehist(
    ax,
    lambda evt: base_mask(evt) & evt.geds.psd.single_temp.is_single_site,
    color="tab:green",
    label="... geds.psd.single_temp.is_single_site",
)
_plot_ehist(
    ax,
    lambda evt: base_mask(evt)
    & evt.geds.psd.single_temp.is_single_site
    & ~evt.coincident.spms,
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
h_spms_mult_sim = hist.new.IntCategory(range(60), name="spms multiplicity").Double()
if has_rc:
    h_spms_mult_rc = hist.new.IntCategory(range(60), name="spms multiplicity").Double()
field_mask = ["spms/multiplicity"]
if has_rc:
    field_mask.append("spms/rc_energy")
it = LH5Iterator(cvt_file, "evt", buffer_len=BUFFER_LEN, field_mask=field_mask)
for evt_chunk in it:
    spms = evt_chunk.view_as("ak").spms
    h_spms_mult_sim.fill(spms.multiplicity)
    if has_rc:
        h_spms_mult_rc.fill(ak.sum(ak.any(spms.rc_energy > 0, axis=-1), axis=-1))
plot.plot_hist(h_spms_mult_sim, ax, label="simulated")
if has_rc:
    plot.plot_hist(h_spms_mult_rc, ax, label="random coincidences")
ax.set_yscale("log")
if has_rc:
    ax.legend()

ax = fig.add_subplot(gs_pe[0, 0])
h_npe_sim = hist.new.Reg(200, 0, 150, name="light per event (photoelectrons)").Double()
if has_rc:
    h_npe_rc = hist.new.Reg(
        200, 0, 150, name="light per event (photoelectrons)"
    ).Double()
field_mask = ["spms/energy_sum"]
if has_rc:
    field_mask.append("spms/rc_energy")
it = LH5Iterator(cvt_file, "evt", buffer_len=BUFFER_LEN, field_mask=field_mask)
for evt_chunk in it:
    spms = evt_chunk.view_as("ak").spms
    h_npe_sim.fill(spms.energy_sum)
    if has_rc:
        h_npe_rc.fill(ak.sum(ak.sum(spms.rc_energy, axis=-1), axis=-1))
plot.plot_hist(h_npe_sim, ax, flow="hint", label="simulated")
if has_rc:
    plot.plot_hist(h_npe_rc, ax, flow="hint", label="random coincidences")
ax.set_ylabel("counts")
ax.set_yscale("log")
ax.legend()

ax = fig.add_subplot(gs_pe[0, 1])
h_time_sim = hist.new.Reg(
    375, -1000, 5000, name="photoelectron $t - t_0$ (ns)"
).Double()
if has_rc:
    h_time_rc = hist.new.Reg(
        375, -1000, 5000, name="photoelectron $t - t_0$ (ns)"
    ).Double()
field_mask = ["spms/time", "trigger/timestamp"]
if has_rc:
    field_mask.append("spms/rc_time")
it = LH5Iterator(cvt_file, "evt", buffer_len=BUFFER_LEN, field_mask=field_mask)
for evt_chunk in it:
    evt = evt_chunk.view_as("ak")
    t0, spms_time = ak.broadcast_arrays(evt.trigger.timestamp, evt.spms.time)
    dt = spms_time - t0
    h_time_sim.fill(ak.flatten(dt, axis=None))
    if has_rc:
        t0_rc, rc_time = ak.broadcast_arrays(evt.trigger.timestamp, evt.spms.rc_time)
        h_time_rc.fill(ak.flatten(rc_time - t0_rc, axis=None))
plot.plot_hist(h_time_sim, ax, flow="none", label="simulated")
if has_rc:
    plot.plot_hist(h_time_rc, ax, flow="none", label="random coincidences")
ax.set_ylabel("counts / 16 ns")
ax.set_yscale("log")
ax.legend()

fig.suptitle(f"{simid}: evt tier")

plot.decorate(fig)
fig.savefig(output_pdf)
