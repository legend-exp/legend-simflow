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


import awkward as ak
import dbetto.utils
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lgdo
import numpy as np
import pyg4ometry
import pygeomhpges
import pygeomtools
import reboost.hpge.psd
import reboost.hpge.surface
import reboost.hpge.utils
import reboost.math.functions
import reboost.spms
from lgdo import lh5
from lgdo.lh5 import LH5Iterator

from legendsimflow import metadata as mutils
from legendsimflow import nersc, patterns, utils
from legendsimflow import reboost as reboost_utils

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

stp_file = args.input.stp_file
jobid = args.wildcards.jobid
hit_file = args.output[0]
gdml_file = args.input.geom
log_file = args.log[0]
metadata = args.config.metadata
hpge_dtmap_files = args.input.hpge_dtmaps
hpge_currmods_files = args.input.hpge_currmods
simstat_part_file = args.input.simstat_part_file
buffer_len = args.params.buffer_len

# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)

# load the geometry and retrieve registered sensitive volumes
geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
sensvols = pygeomtools.detectors.get_all_sensvols(geom)

partitions = dbetto.utils.load_dict(simstat_part_file)[f"job_{jobid}"]


# loop over the partitions for this file
for runid, evt_idx_range in partitions.items():
    msg = f"processing partition corresponding to {runid}, event range {evt_idx_range}"
    log.info(msg)

    # load in parameters of the HPGe current signal model
    currmod_pars_file = patterns.output_currmod_merged_filename(
        snakemake.config,  # noqa: F821
        runid=runid,
    )
    currmod_pars_all = dbetto.utils.load_dict(currmod_pars_file)

    # loop over the sensitive volumes registered in the geometry
    for det_name, geom_meta in sensvols.items():
        msg = f"looking for data from sensitive volume {det_name} (uid={geom_meta.uid})..."
        log.debug(msg)

        # get the usability
        usability = mutils.usability(metadata, det_name, runid=runid, default="on")

        if f"stp/{det_name}" not in lh5.ls(stp_file, "*/*"):
            msg = (
                f"detector {det_name} not found in {stp_file}. "
                "possibly because it was not read-out or there were no hits recorded"
            )
            log.warning(msg)
            continue

        i_start, n_entries = reboost_utils.get_remage_hit_range(
            stp_file, det_name, geom_meta.uid, evt_idx_range
        )

        # initialize the stp file iterator
        # NOTE: if the entry list is empty, there will be no processing but an
        # empty output table will be nonetheless created. this is important for
        # the buil_tcm() step at the end
        iterator = LH5Iterator(
            stp_file,
            f"stp/{det_name}",
            i_start=i_start,
            n_entries=n_entries,
            buffer_len=buffer_len,
        )

        # only process the HPGe output
        if geom_meta.detector_type != "germanium":
            continue

        msg = f"processing the {det_name} output table..."
        log.info(msg)

        log.debug("creating an pygeomhpges.HPGe object")
        pyobj = pygeomhpges.make_hpge(
            geom_meta.metadata, registry=None, allow_cylindrical_asymmetry=False
        )

        fccd = mutils.get_sanitized_fccd(metadata, det_name)
        det_loc = geom.physicalVolumeDict[det_name].position

        # NOTE: we don't use the script arg but we use the (known) file patterns. more robust
        dt_map = reboost_utils.load_hpge_dtmaps(snakemake.config, det_name, runid)  # noqa: F821

        # load parameters of the current model
        pars = currmod_pars_all.get(det_name, None)
        currmod_pars = (
            pars.get("current_pulse_shape", None) if pars is not None else None
        )

        # iterate over input data
        for lgdo_chunk in iterator:
            chunk = lgdo_chunk.view_as("ak")

            _distance_to_nplus = reboost.hpge.surface.distance_to_surface(
                chunk.xloc * 1000,  # mm
                chunk.yloc * 1000,  # mm
                chunk.zloc * 1000,  # mm
                pyobj,
                det_loc.eval(),
                distances_precompute=chunk.dist_to_surf * 1000,
                precompute_cutoff=(fccd + 1),
                surface_type="nplus",
            )

            _activeness = reboost.math.functions.piecewise_linear_activeness(
                _distance_to_nplus,
                fccd=fccd,
                dlf=0.5,
            )

            edep_active = chunk.edep * _activeness
            energy = ak.sum(edep_active, axis=-1)

            # PSD: if the drift time map is none, it means that we don't
            # have the detector model to simulate PSD in a more advanced
            # way

            # default to NaN
            dt_heuristic = np.full(len(chunk), np.nan)
            aoe = np.full(len(chunk), np.nan)

            if dt_map is not None:
                _drift_time = reboost_utils.hpge_corrected_drift_time(
                    chunk, dt_map, det_loc
                )
                dt_heuristic = reboost.hpge.psd.drift_time_heuristic(
                    _drift_time, chunk.edep
                )
                _a_max = reboost_utils.hpge_max_current_cal(
                    edep_active, _drift_time, currmod_pars
                )
                aoe = _a_max / energy

            out_table = reboost_utils.make_output_chunk(lgdo_chunk)
            out_table.add_field("energy", lgdo.Array(energy, attrs={"units": "keV"}))
            out_table.add_field("drift_time_heuristic", lgdo.Array(dt_heuristic))
            out_table.add_field("aoe", lgdo.Array(aoe))

            # add strings
            utils.add_field_string("runid", out_table, runid)
            utils.add_field_string("usability", out_table, usability)

            reboost_utils.write_chunk(
                out_table,
                f"/hit/{det_name}",
                hit_file,
                geom_meta.uid,
            )

log.debug("building the TCM")
reboost_utils.build_tcm(hit_file)
