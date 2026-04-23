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

import argparse
from pathlib import Path

import awkward as ak
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import numpy as np
import pyg4ometry
import pygeomtools
import reboost.hpge.psd
import reboost.hpge.surface
import reboost.hpge.utils
import reboost.math.functions
import reboost.spms
from dbetto import AttrsDict
from dbetto.utils import load_dict
from lgdo import Array, VectorOfVectors, lh5
from lgdo.lh5 import LH5Iterator
from reboost.optmap.convolve import OptmapForConvolve
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import metadata as mutils
from legendsimflow import nersc, utils
from legendsimflow import reboost as reboost_utils
from legendsimflow.metadata import get_tier_settings
from legendsimflow.profile import make_profiler
from legendsimflow.scripts import log_script_invocation
from legendsimflow.tcm import build_tcm


@snakemake_compatible(
    mapping={
        "stp_file": "input.stp_file",
        "optmap_lar": "input.optmap_lar",
        "geom_file": "input.geom",
        "simstat_part_file": "input.simstat_part_file",
        "detector_usabilities_file": "input.detector_usabilities[0]",
        "jobid": "wildcards.jobid",
        "opt_file": "output[0]",
        "log_file": "log[0]",
        "optmap_per_sipm": "params.optmap_per_sipm",
        "scintillator_volume_name": "params.scintillator_volume_name",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(description="Build the opt tier.")
    parser.add_argument("--stp-file", required=True, help="input stp tier file")
    parser.add_argument("--optmap-lar", required=True, help="LAr optical map file")
    parser.add_argument("--geom-file", required=True, help="GDML geometry file")
    parser.add_argument(
        "--simstat-part-file",
        required=True,
        help="simulation statistics partition file",
    )
    parser.add_argument(
        "--detector-usabilities-file",
        required=True,
        help="detector usabilities YAML file",
    )
    parser.add_argument("--jobid", required=True, help="job ID wildcard")
    parser.add_argument("--opt-file", required=True, help="output opt tier file")
    parser.add_argument("--log-file", default=None, help="log file")
    parser.add_argument(
        "--optmap-per-sipm",
        action="store_true",
        default=False,
        help="sample photoelectrons per SiPM (default: all SiPMs together)",
    )
    parser.add_argument(
        "--scintillator-volume-name",
        required=True,
        help="name of the scintillator sensitive volume in the geometry",
    )
    parser.add_argument(
        "--simflow-config",
        "--config",
        dest="simflow_config",
        required=True,
        help="simflow config YAML path",
    )
    args = parser.parse_args()

    config = utils.init_simflow_context(args.simflow_config, workflow=None).config

    stp_file = nersc.dvs_ro(config, args.stp_file)
    jobid = args.jobid
    opt_file = args.opt_file
    optmap_lar = nersc.dvs_ro(config, args.optmap_lar)
    gdml_file = nersc.dvs_ro(config, args.geom_file)
    log_file = args.log_file
    metadata = config.metadata
    optmap_per_sipm = args.optmap_per_sipm
    scintillator_volume_name = args.scintillator_volume_name
    simstat_part_file = nersc.dvs_ro(config, args.simstat_part_file)
    usabilities = AttrsDict(
        load_dict(nersc.dvs_ro(config, args.detector_usabilities_file))
    )

    opt_file, move2cfs = nersc.make_on_scratch(config, opt_file)

    tier_opt_settings = get_tier_settings(config, "opt")
    optmap_scaling_factor = tier_opt_settings.optmap_scaling_factor
    photoelectron_resolution_sigma = tier_opt_settings.photoelectron_resolution_sigma
    time_resolution_in_ns = tier_opt_settings.time_resolution_in_ns
    max_pes_per_hit = (
        tier_opt_settings.max_pes_per_hit_per_sipm
        if optmap_per_sipm
        else tier_opt_settings.max_pes_per_hit_combined
    )
    buffer_len = tier_opt_settings.buffer_len

    # setup logging
    log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "tier-opt", parser, args)
    perf_block, print_perf, print_perf_last = make_profiler()

    # load the geometry and retrieve registered sensitive volume tables
    geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
    sens_tables = pygeomtools.detectors.get_all_senstables(geom)

    def process_sipm(
        iterator: LH5Iterator,
        optmap_lar: str | Path | OptmapForConvolve,
        sipm: str,
        sipm_uid: int,
        out_file: str | Path,
        runid: str,
        usability: str,
    ) -> None:
        with perf_block("load_optmap()"):
            if not isinstance(optmap_lar, OptmapForConvolve):
                optmap_lar = reboost.spms.pe.load_optmap(optmap_lar, sipm)

        for lgdo_chunk in iterator:
            chunk = lgdo_chunk.view_as("ak")

            with perf_block("emitted_scintillation_photons()"):
                scint_ph = reboost.spms.pe.emitted_scintillation_photons(
                    chunk.edep, chunk.particle, "lar"
                )

            with perf_block("number_of_detected_photoelectrons()"):
                _output = reboost.spms.pe.number_of_detected_photoelectrons(
                    chunk.xloc,
                    chunk.yloc,
                    chunk.zloc,
                    scint_ph,
                    optmap_lar,
                    sipm,
                    map_scaling=optmap_scaling_factor,
                    max_pes_per_hit=max_pes_per_hit,
                )
            if max_pes_per_hit > 0:
                nr_pe, is_saturated = _output
            else:
                nr_pe = _output
                is_saturated = np.full(len(chunk), fill_value=False, dtype=np.bool_)

            with perf_block("photoelectron_times()"):
                pe_times_micro = reboost.spms.pe.photoelectron_times(
                    nr_pe, chunk.particle, chunk.time, "lar"
                )

                # the photoelectron_times() processor does not guarantee time
                # ordering
                pe_times_micro = ak.sort(pe_times_micro, axis=-1)

            with perf_block("photoelectron_resolution()"):
                pe_amps_micro = reboost_utils.smear_photoelectrons(
                    pe_times_micro, photoelectron_resolution_sigma
                )

            if time_resolution_in_ns > 0:
                with perf_block("cluster_photoelectrons()"):
                    pe_times, pe_amps = reboost_utils.cluster_photoelectrons(
                        pe_times_micro,
                        pe_amps_micro,
                        time_resolution_in_ns,
                    )
            else:
                pe_times = pe_times_micro
                pe_amps = pe_amps_micro

            with perf_block("write_chunk()"):
                out_table = reboost_utils.make_output_chunk(lgdo_chunk)

                out_table.add_field(
                    "time",
                    VectorOfVectors(
                        ak.values_astype(pe_times, np.float32), attrs={"units": "ns"}
                    ),
                )
                out_table.add_field(
                    "energy", VectorOfVectors(ak.values_astype(pe_amps, np.float32))
                )
                out_table.add_field("is_saturated", Array(is_saturated))

                _, period, run, _ = mutils.parse_runid(runid)
                field_vals = [period, run, mutils.encode_usability(usability)]
                for i, field in enumerate(["period", "run", "usability"]):
                    out_table.add_field(
                        field,
                        Array(np.full(shape=len(chunk), fill_value=field_vals[i])),
                    )

                reboost_utils.write_chunk(
                    out_table,
                    "/hit/" + ("spms" if sipm == "all" else sipm),
                    out_file,
                    sipm_uid,
                )

    partitions = load_dict(simstat_part_file)[f"job_{jobid}"]

    # load TCM, to be used to chunk the event statistics according to the run partitioning
    msg = "loading TCM"
    log.debug(msg)
    tcm = lh5.read_as("tcm", stp_file, library="ak")

    # pre load optical map for a little speed up
    if not optmap_per_sipm:
        optmap_lar = reboost.spms.pe.load_optmap(optmap_lar, "all")

    # loop over the partitions for this file
    for runid_idx, (runid, evt_idx_range) in enumerate(partitions.items()):
        msg = (
            f"processing partition corresponding to {runid} "
            f"[{runid_idx + 1}/{len(partitions)}], event range {evt_idx_range}"
        )
        log.info(msg)

        # loop over the sensitive volume tables registered in the geometry
        for det_name, geom_meta in sens_tables.items():
            # process the scintillator output
            if not (
                geom_meta.detector_type == "scintillator"
                and det_name == scintillator_volume_name
            ):
                continue

            msg = f"looking for data from sensitive volume {det_name} table (uid={geom_meta.uid})..."
            log.debug(msg)

            if f"stp/{det_name}" not in lh5.ls(stp_file, "*/*"):
                msg = (
                    f"detector {det_name} not found in {stp_file}. "
                    "possibly because it was not read-out or there were no hits recorded"
                )
                log.warning(msg)
                continue

            log.info("processing the 'lar' scintillator table...")

            msg = "looking for indices of hit table rows to read..."
            log.debug(msg)
            i_start, n_entries = reboost_utils.get_remage_hit_range(
                tcm, det_name, geom_meta.uid, evt_idx_range
            )

            def _make_iterator(det_name=det_name, i_start=i_start, n_entries=n_entries):
                return LH5Iterator(
                    stp_file,
                    f"stp/{det_name}",
                    i_start=i_start,
                    n_entries=n_entries,
                    buffer_len=buffer_len,
                )

            if optmap_per_sipm:
                for sipm in sorted(reboost_utils.get_senstables(geom, "optical")):
                    sipm_uid = sens_tables[sipm].uid

                    # get the usability
                    det_info = usabilities[runid].get(sipm, None)
                    if det_info is None:
                        msg = f"usability not found for {sipm} in {runid}, defaulting to on"
                        log.warning(msg)
                        usability = "on"
                    else:
                        usability = det_info.usability

                    msg = f"applying optical map for SiPM {sipm}"
                    log.debug(msg)

                    process_sipm(
                        _make_iterator(),
                        optmap_lar,
                        sipm,
                        sipm_uid,
                        opt_file,
                        runid,
                        usability,
                    )

                    print_perf_last()
            else:
                log.debug("applying sum optical map")

                process_sipm(
                    _make_iterator(),
                    optmap_lar,
                    "all",
                    geom_meta.uid,
                    opt_file,
                    runid,
                    "on",
                )

    log.debug("building the TCM")
    build_tcm(opt_file, opt_file)

    with perf_block("move_to_cfs()"):
        move2cfs()

    print_perf()


if __name__ == "__main__":
    main()
