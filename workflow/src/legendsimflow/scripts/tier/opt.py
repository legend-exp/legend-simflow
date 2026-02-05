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

from pathlib import Path

import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import pyg4ometry
import pygeomtools
import reboost.hpge.psd
import reboost.hpge.surface
import reboost.hpge.utils
import reboost.math.functions
import reboost.spms
from lgdo import VectorOfVectors, lh5
from lgdo.lh5 import LH5Iterator

from legendsimflow import nersc
from legendsimflow import reboost as reboost_utils
from legendsimflow.profile import make_profiler

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

stp_file = args.input.stp_file
hit_file = args.output[0]
optmap_lar_file = args.input.optmap_lar
gdml_file = args.input.geom
log_file = args.log[0]
metadata = args.config.metadata
optmap_per_sipm = args.params.optmap_per_sipm
scintillator_volume_name = args.params.scintillator_volume_name

# for some sims like Th228 loading a 100MB chunk of the TCM can result in a lot
# of photons, i.e. high memory usage
BUFFER_LEN = "50*MB"
MAP_SCALING = 0.1
DEFAULT_PHOTOELECTRON_RES = 0.2  # FWHM
TIME_RESOLUTION_NS = 16

# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
perf_block, print_perf = make_profiler()

# load the geometry and retrieve registered sensitive volume tables
geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
sens_tables = pygeomtools.detectors.get_all_senstables(geom)


def process_sipm(
    iterator: LH5Iterator,
    optmap_lar_file: str | Path,
    sipm: str,
    sipm_uid: int,
    out_file: str | Path,
) -> None:
    with perf_block("load_optmap()"):
        optmap = reboost.spms.pe.load_optmap(optmap_lar_file, sipm)

    for lgdo_chunk in iterator:
        chunk = lgdo_chunk.view_as("ak")

        with perf_block("emitted_scintillation_photons()"):
            scint_ph = reboost.spms.pe.emitted_scintillation_photons(
                chunk.edep, chunk.particle, "lar"
            )

        with perf_block("number_of_detected_photoelectrons()"):
            nr_pe = reboost.spms.pe.number_of_detected_photoelectrons(
                chunk.xloc,
                chunk.yloc,
                chunk.zloc,
                scint_ph,
                optmap,
                sipm,
                map_scaling=MAP_SCALING,
            )

        with perf_block("photoelectron_times()"):
            pe_times_micro = reboost.spms.pe.photoelectron_times(
                nr_pe, chunk.particle, chunk.time, "lar"
            ).view_as("ak")

        with perf_block("photoelectron_resolution()"):
            pe_amps_micro = reboost_utils.smear_photoelectrons(
                pe_times_micro, DEFAULT_PHOTOELECTRON_RES
            )

        with perf_block("cluster_photoelectrons()"):
            pe_times, pe_amps = reboost_utils.cluster_photoelectrons(
                pe_times_micro,
                pe_amps_micro,
                TIME_RESOLUTION_NS,
            )

        with perf_block("write_chunk()"):
            out_table = reboost_utils.make_output_chunk(lgdo_chunk)
            out_table.add_field(
                "time", VectorOfVectors(pe_times, attrs={"units": "ns"})
            )
            out_table.add_field("energy", VectorOfVectors(pe_amps))
            reboost_utils.write_chunk(
                out_table,
                "/hit/" + ("spms" if sipm == "all" else sipm),
                out_file,
                sipm_uid,
            )


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

    # process the scintillator output
    if (
        geom_meta.detector_type == "scintillator"
        and det_name == scintillator_volume_name
    ):
        log.info("processing the 'lar' scintillator table...")

        if optmap_per_sipm:
            for sipm in reboost_utils.get_senstables(geom, "optical"):
                sipm_uid = sens_tables[sipm].uid

                msg = f"applying optical map for SiPM {sipm}"
                log.debug(msg)

                iterator = LH5Iterator(
                    stp_file,
                    f"stp/{det_name}",
                    buffer_len=BUFFER_LEN,
                )

                process_sipm(iterator, optmap_lar_file, sipm, sipm_uid, hit_file)

        else:
            log.debug("applying sum optical map")

            iterator = LH5Iterator(
                stp_file,
                f"stp/{det_name}",
                buffer_len=BUFFER_LEN,
            )

            process_sipm(iterator, optmap_lar_file, "all", geom_meta.uid, hit_file)


log.debug("building the TCM")
reboost_utils.build_tcm(hit_file, hit_file)

print_perf()
