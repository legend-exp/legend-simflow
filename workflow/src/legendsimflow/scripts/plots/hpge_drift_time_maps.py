# ruff: noqa: I002

# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>,
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

import matplotlib.pyplot as plt
import numpy as np
from lgdo import lh5

from legendsimflow import nersc
from legendsimflow.plot import decorate

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

dtmap_file = args.input[0]
output_pdf = args.output[0]
hpge = args.wildcards.hpge_detector

dtmap = lh5.read(hpge, dtmap_file)

r = dtmap.r.view_as("np")
zcoord = dtmap.z.view_as("np")


def symmetrize(a):
    a = a.T
    return np.concatenate((np.fliplr(a), a), axis=1)


img_045 = symmetrize(dtmap.drift_time_045_deg.view_as("np"))
img_000 = symmetrize(dtmap.drift_time_000_deg.view_as("np"))

# ratio (mask invalid divisions)
ratio = np.divide(
    img_000,
    img_045,
    out=np.full_like(img_000, np.nan),
    where=img_045 > 0,
)

extent = (-r.max(), r.max(), zcoord.min(), zcoord.max())

vmin = np.nanmin([img_045, img_000])
vmax = np.nanmax([img_045, img_000])

fig, axes = plt.subplots(
    ncols=3,
    figsize=(12, 4),
    sharey=True,
)


def plot(ax, img, title, *, cmap, vmin=None, vmax=None):
    im = ax.imshow(
        img,
        origin="lower",
        extent=extent,
        aspect="equal",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_xlabel("r (m)")
    ax.set_title(f"{hpge} · {title}")
    return im


im0 = plot(axes[0], img_045, "<110>", cmap="viridis", vmin=vmin, vmax=vmax)
axes[0].set_ylabel("z (m)")

im1 = plot(axes[1], img_000, "<001>", cmap="viridis", vmin=vmin, vmax=vmax)

im2 = plot(axes[2], ratio, "<001> / <110>", cmap="coolwarm")

fig.colorbar(im1, ax=axes[:2], label="drift time (ns)")
fig.colorbar(im2, ax=axes[2], label="ratio")

decorate(fig)

plt.savefig(output_pdf)
