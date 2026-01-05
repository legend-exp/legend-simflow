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

remage_file = args.input[0]
output_pdf = args.output[0]

vtx = lh5.read_as("/vtx", remage_file, n_rows=10_000, library="ak")

x = vtx.xloc.to_numpy()
y = vtx.yloc.to_numpy()
z = vtx.zloc.to_numpy()

fig = plt.figure(figsize=(16, 6))
gs = fig.add_gridspec(1, 3, wspace=0.1)

ax_xy = fig.add_subplot(gs[0, 0])
ax_3d = fig.add_subplot(gs[0, 1], projection="3d")
ax_xz = fig.add_subplot(gs[0, 2])

# set consistent limits using a single scale across x,y,z
mins = np.array([x.min(), y.min(), z.min()], dtype=float)
maxs = np.array([x.max(), y.max(), z.max()], dtype=float)
centers = 0.5 * (mins + maxs)
half_range = 0.5 * (maxs - mins).max()
half_range *= 1.1

xlim = (centers[0] - half_range, centers[0] + half_range)
ylim = (centers[1] - half_range, centers[1] + half_range)
zlim = (centers[2] - half_range, centers[2] + half_range)

# apply limits everywhere
ax_xy.set_xlim(xlim)
ax_xy.set_ylim(ylim)

ax_xz.set_xlim(xlim)
ax_xz.set_ylim(zlim)

ax_3d.set_xlim(xlim)
ax_3d.set_ylim(ylim)
ax_3d.set_zlim(zlim)

# enforce equal scaling
ax_xy.set_aspect("equal", adjustable="box")
ax_xz.set_aspect("equal", adjustable="box")
ax_3d.set_box_aspect((1, 1, 1))  # equal x:y:z in display space

# plots
ax_xy.scatter(x, y, s=1)
ax_xy.set_xlabel("x [m]")
ax_xy.set_ylabel("y [m]")
ax_xy.set_title("xy projection")

ax_3d.scatter(x, y, z, s=1)
ax_3d.set_xlabel("x [m]")
ax_3d.set_ylabel("y [m]")
ax_3d.set_zlabel("z [m]")
ax_3d.set_title("3D")

ax_xz.scatter(x, z, s=1)
ax_xz.set_xlabel("x [m]")
ax_xz.set_ylabel("z [m]")
ax_xz.set_title("xz projection")

for ax in (ax_xy, ax_xz):
    ax.grid(True)

decorate(fig)
plt.savefig(output_pdf)
