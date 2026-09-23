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

import lh5
import matplotlib.pyplot as plt
import pyg4ometry
import pygeomhpges
from matplotlib.backends.backend_pdf import PdfPages
from reboost.hpge import plot_rz_maps

from legendsimflow import nersc
from legendsimflow.plot import decorate

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

dtmap_file = args.input[0]
output_pdf = args.output[0]

reg = pyg4ometry.geant4.Registry()
natge = pygeomhpges.materials.make_natural_germanium(registry=reg)


def fig(hpge):
    # HPGe profile
    pyobj = pygeomhpges.make_hpge(
        args.config.metadata.hardware.detectors.germanium.diodes[hpge],
        registry=reg,
        material=natge,
        allow_cylindrical_asymmetry=False,
    )

    dtmap = lh5.read(hpge, dtmap_file)

    # grid in m in the file, in mm in pygeomhpges; maps keyed by crystal axis azimuth
    fig, _ = plot_rz_maps(
        {
            0: dtmap.drift_time_000_deg.view_as("np"),
            45: dtmap.drift_time_045_deg.view_as("np"),
        },
        dtmap.r.view_as("np") * 1000,
        dtmap.z.view_as("np") * 1000,
        hpge=pyobj,
        label="drift time [ns]",
        title=hpge,
    )

    return fig


# prepare a pdf file with a plot per page
tables = sorted(lh5.ls(dtmap_file))
fig_builders = [lambda t=t: fig(t) for t in tables]

with PdfPages(output_pdf) as pdf:
    for make_fig in fig_builders:
        fig = make_fig()
        decorate(fig)
        pdf.savefig(fig)
        plt.close(fig)
