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
import pygeomhpges
import pygeomhpges.draw
from lgdo import lh5
from matplotlib.backends.backend_pdf import PdfPages

from legendsimflow import nersc
from legendsimflow.plot import decorate

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

dtmap_file = args.input[0]
output_pdf = args.output[0]


def symmetrize(a):
    a = a.T
    return np.concatenate((np.fliplr(a), a), axis=1)


def save_page(pdf, make_fig):
    fig = make_fig()
    decorate(fig)
    pdf.savefig(fig)
    plt.close(fig)


def fig(hpge):
    # HPGe profile
    pyobj = pygeomhpges.make_hpge(
        args.config.metadata.hardware.detectors.germanium.diodes[hpge], registry=None
    )

    dtmap = lh5.read(hpge, dtmap_file)

    r = dtmap.r.view_as("np") * 1000  # convert to mm, units in pygeomhpges
    zcoord = dtmap.z.view_as("np") * 1000

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
        pygeomhpges.draw.plot_profile(
            pyobj,
            axes=ax,
            marker=None,
            linewidth=1,
            color="black",
        )

        xmin, xmax = -r.max(), r.max()
        ymin, ymax = zcoord.min(), zcoord.max()

        f = 0.04
        ax.set_xlim(xmin - f * (xmax - xmin), xmax + f * (xmax - xmin))
        ax.set_ylim(ymin - f * (ymax - ymin), ymax + f * (ymax - ymin))

        ax.set_xlabel("r (mm)")
        ax.set_title(f"{hpge} · {title}")

        return im

    plot(axes[0], img_045, "<110>", cmap="viridis", vmin=vmin, vmax=vmax)
    axes[0].set_ylabel("z (mm)")

    im1 = plot(axes[1], img_000, "<001>", cmap="viridis", vmin=vmin, vmax=vmax)

    im2 = plot(axes[2], ratio, "<001> / <110>", cmap="coolwarm")

    fig.colorbar(im1, ax=axes[:2], label="drift time (ns)")
    fig.colorbar(im2, ax=axes[2], label="ratio")

    return fig


# prepare a pdf file with a plot per page
tables = sorted(lh5.ls(dtmap_file))
fig_builders = [lambda t=t: fig(t) for t in tables]

with PdfPages(output_pdf) as pdf:
    for make_fig in fig_builders:
        save_page(pdf, make_fig)
