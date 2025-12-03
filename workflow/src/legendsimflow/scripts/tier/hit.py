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
import legendhpges
import lgdo
import numpy as np
import pyg4ometry
import pygama.evt
import pygeomtools
import reboost.hpge.psd
import reboost.hpge.surface
import reboost.hpge.utils
import reboost.math.functions
import reboost.spms
from dbetto import AttrsDict
from legendmeta.police import validate_dict_schema
from lgdo import lh5
from lgdo.lh5 import LH5Iterator

from legendsimflow import reboost as reboost_utils

args = snakemake  # noqa: F821

stp_file = args.input.stp_file
jobid = args.wildcards.jobid
hit_file = args.output[0]
gdml_file = args.input.geom
log_file = args.log[0]
metadata = args.config.metadata
hpge_dtmap_files = args.input.hpge_dtmaps
simstat_part_file = args.input.simstat_part_file
buffer_len = args.params.buffer_len


# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)

# load the geometry and retrieve registered sensitive volumes
geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
sensvols = pygeomtools.detectors.get_all_sensvols(geom)

# load TCM, to be used to chunk the event statistics according to the run partitioning
tcm = lh5.read_as("tcm", stp_file, library="ak")
partitions = dbetto.utils.load_dict(simstat_part_file)[f"job_{jobid}"]


# loop over the partitions for this file
for runid, tcm_idx_range in partitions.items():
    msg = f"processing partition corresponding to {runid}, event range {tcm_idx_range}"
    log.info(msg)

    # loop over the sensitive volumes registered in the geometry
    for det_name, geom_meta in sensvols.items():
        msg = f"looking for data from sensitive volume {det_name} (uid={geom_meta.uid})..."
        log.debug(msg)

        if f"stp/{det_name}" not in lh5.ls(stp_file, "*/*"):
            msg = (
                f"detector {det_name} not found in {stp_file}. "
                "possibly because it was not read-out or there were no hits recorded"
            )
            log.warning(msg)
            continue

        # ask the TCM which rows we should read from the hit table
        tcm_part = tcm[tcm_idx_range[0] : tcm_idx_range[1]]
        entry_list = ak.flatten(
            tcm_part[tcm_part.table_key == geom_meta.uid].row_in_table
        ).to_list()

        if len(entry_list) > 0:
            assert list(range(entry_list[0], entry_list[-1] + 1)) == entry_list

            msg = (
                f"hits with indices in [{entry_list[0]}, {entry_list[-1]}] "
                "recorded in the events belonging to this partition"
            )
            log.debug(msg)

            i_start = entry_list[0]
            n_entries = entry_list[-1] - entry_list[0]

        else:
            msg = f"no hits recorded in {det_name} in the events belonging to this partition"
            log.warning(msg)

            i_start = 0
            n_entries = None

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

        # process the HPGe output
        if geom_meta.detector_type == "germanium":
            msg = f"processing the {det_name} output table..."
            log.info(msg)

            log.debug("creating an legendhpges.HPGe object")
            pyobj = legendhpges.make_hpge(
                geom_meta.metadata, registry=None, allow_cylindrical_asymmetry=False
            )

            det_meta = metadata.hardware.detectors.germanium.diodes[det_name]

            has_fccd_meta = validate_dict_schema(
                det_meta.characterization,
                {"combined_0vbb_analysis": {"fccd_in_mm": 0}},
                greedy=False,
                verbose=False,
            )

            if not has_fccd_meta:
                msg = f"{det_name} metadata does not seem to contain usable FCCD data, setting to 1 mm"
                log.warning(msg)
                fccd = 1
            else:
                fccd = det_meta.characterization.combined_0vbb_analysis.fccd_in_mm

            det_loc = geom.physicalVolumeDict[det_name].position

            dt_map = reboost_utils.load_hpge_dtmaps(snakemake.config, det_name, runid)  # noqa: F821

            # iterate over input data
            for lgdo_chunk in iterator:
                chunk = lgdo_chunk.view_as("ak")
                _distance_to_surf = AttrsDict()

                for surf in ("nplus", "pplus", "passive"):
                    _distance_to_surf[surf] = reboost.hpge.surface.distance_to_surface(
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
                    _distance_to_surf.nplus,
                    fccd=fccd,
                    dlf=0.2,
                )

                energy = ak.sum(chunk.edep * _activeness, axis=-1)

                if dt_map is not None:
                    _phi = np.arctan2(
                        chunk.yloc * 1000 - det_loc.eval()[1],
                        chunk.xloc * 1000 - det_loc.eval()[0],
                    )

                    _drift_time = {}
                    for angle, _map in dt_map.items():
                        _drift_time[angle] = reboost.hpge.psd.drift_time(
                            chunk.xloc,
                            chunk.yloc,
                            chunk.zloc,
                            _map,
                            coord_offset=det_loc,
                        ).view_as("ak")

                    _drift_time_corr = (
                        _drift_time["045"]
                        + (_drift_time["000"] - _drift_time["045"])
                        * (1 - np.cos(4 * _phi))
                        / 2
                    )

                    dt_heuristic = reboost.hpge.psd.drift_time_heuristic(
                        _drift_time_corr, chunk.edep
                    )
                else:
                    dt_heuristic = np.full(len(chunk), np.nan)

                out_table = reboost_utils.make_output_chunk(lgdo_chunk)
                out_table.add_field(
                    "energy", lgdo.Array(energy, attrs={"units": "keV"})
                )
                out_table.add_field("drift_time_heuristic", lgdo.Array(dt_heuristic))

                reboost_utils.write_chunk(
                    out_table,
                    f"/hit/{det_name}",
                    hit_file,
                    geom_meta.uid,
                    runid,
                )

# build the TCM
# use tables keyed by UID in the __by_uid__ group.  in this way, the
# TCM will index tables by UID.  the coincidence criterium is based
# on Geant4 event identifier and time of the hits
# NOTE: uses the same time window as in build_hit() reshaping
log.debug("building the TCM")
pygama.evt.build_tcm(
    [(hit_file, r"hit/__by_uid__/*")],  # input_tables
    ["evtid", "t0"],  # coin_cols
    hash_func=r"(?<=hit/__by_uid__/det)\d+",
    coin_windows=[0, 10_000],
    out_file=hit_file,
    wo_mode="write_safe",
)
