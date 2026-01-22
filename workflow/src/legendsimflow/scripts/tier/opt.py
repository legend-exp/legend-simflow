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

import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import pyg4ometry
import pygeomtools
import reboost.hpge.psd
import reboost.hpge.surface
import reboost.hpge.utils
import reboost.math.functions
import reboost.spms
from lgdo import lh5
from lgdo.lh5 import LH5Iterator

from legendsimflow import nersc
from legendsimflow import reboost as reboost_utils

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

stp_file = args.input.stp_file
hit_file = args.output[0]
optmap_lar_file = args.input.optmap_lar
gdml_file = args.input.geom
log_file = args.log[0]
metadata = args.config.metadata
optmap_per_sipm = args.params.optmap_per_sipm
buffer_len = "100*MB"
scintillator_volume_name = args.params.scintillator_volume_name


# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)

# load the geometry and retrieve registered sensitive volume tables
geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
sens_tables = pygeomtools.detectors.get_all_senstables(geom)

# loop over the sensitive volume tables registered in the geometry
for det_name, geom_meta in sens_tables.items():
    msg = f"looking for data from sensitive volume {det_name} table (uid={geom_meta.uid})..."
    log.debug(msg)

    if f"stp/{det_name}" not in lh5.ls(stp_file, "*/*"):
        msg = (
            f"detector {det_name} not found in {stp_file}. "
            "possibly because it was not read-out or there were no hits recorded"
        )
        log.warning(msg)
        continue

    # initialize the stp file iterator
    # NOTE: if the entry list is empty, there will be no processing but an
    # empty output table will be nonetheless created
    iterator = LH5Iterator(
        stp_file,
        f"stp/{det_name}",
        buffer_len=buffer_len,
    )

    # process the scintillator output
    if (
        geom_meta.detector_type == "scintillator"
        and det_name == scintillator_volume_name
    ):
        log.info("processing the 'lar' scintillator table...")

        # QUESTION/FIXME: what is the right loop order? load a map and process
        # all chunks or load a chunk and loop over the maps?
        if not optmap_per_sipm:
            optmap = reboost.spms.pe.load_optmap(optmap_lar_file, "all")

        for lgdo_chunk in iterator:
            chunk = lgdo_chunk.view_as("ak")

            _scint_ph = reboost.spms.pe.emitted_scintillation_photons(
                chunk.edep, chunk.particle, "lar"
            )

            if optmap_per_sipm:
                for sipm in reboost_utils.get_senstables(geom, "optical"):
                    sipm_uid = sens_tables[sipm].uid

                    msg = f"applying optical map for SiPM {sipm}"
                    log.debug(msg)

                    optmap = reboost.spms.pe.load_optmap(optmap_lar_file, sipm)

                    photoelectrons = reboost.spms.pe.detected_photoelectrons(
                        _scint_ph,
                        chunk.particle,
                        chunk.time,
                        chunk.xloc,
                        chunk.yloc,
                        chunk.zloc,
                        optmap,
                        "lar",
                        sipm,
                        map_scaling=0.1,
                    )

                    out_table = reboost_utils.make_output_chunk(lgdo_chunk)
                    out_table.add_field("time", photoelectrons)
                    reboost_utils.write_chunk(
                        out_table,
                        f"/hit/{sipm}",
                        hit_file,
                        sipm_uid,
                    )
            else:
                log.debug("applying sum optical map")

                photoelectrons = reboost.spms.pe.detected_photoelectrons(
                    _scint_ph,
                    chunk.particle,
                    chunk.time,
                    chunk.xloc,
                    chunk.yloc,
                    chunk.zloc,
                    optmap,
                    "lar",
                    "all",
                    map_scaling=0.1,
                )

                out_table = reboost_utils.make_output_chunk(lgdo_chunk)
                out_table.add_field("time", photoelectrons)
                reboost_utils.write_chunk(
                    out_table,
                    "/hit/spms",
                    hit_file,
                    geom_meta.uid,
                )


log.debug("building the TCM")
reboost_utils.build_tcm(hit_file, hit_file)
