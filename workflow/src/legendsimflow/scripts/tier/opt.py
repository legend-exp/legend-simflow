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

import random
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

from legendsimflow import metadata as mutils
from legendsimflow import nersc, spms_pars, utils
from legendsimflow import reboost as reboost_utils
from legendsimflow.profile import make_profiler
from legendsimflow.tcm import build_tcm

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

stp_file = args.input.stp_file
jobid = args.wildcards.jobid
hit_file = args.output[0]
optmap_lar = args.input.optmap_lar
gdml_file = args.input.geom
log_file = args.log[0]
metadata = args.config.metadata
optmap_per_sipm = args.params.optmap_per_sipm
add_random_coincidences = args.params.add_random_coincidences
scintillator_volume_name = args.params.scintillator_volume_name
simstat_part_file = args.input.simstat_part_file
usabilities = AttrsDict(load_dict(args.input.detector_usabilities[0]))
l200data = args.config.paths.l200data

hit_file, move2cfs = nersc.make_on_scratch(args.config, hit_file)

# for some sims like Th228 loading a 100MB chunk of the TCM can result in a lot
# of photons, i.e. high memory usage
BUFFER_LEN = "10*MB"
MAP_SCALING = 1  # FIXME: guess
DEFAULT_PHOTOELECTRON_RES = 0.3  # FWHM FIXME: guess
TIME_RESOLUTION_NS = 16  # FIXME: guess
MAX_PES_PER_HIT = 5 if optmap_per_sipm else 100

# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
perf_block, print_perf, print_perf_last = make_profiler()

# load the geometry and retrieve registered sensitive volume tables
geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
sens_tables = pygeomtools.detectors.get_all_senstables(geom)


def _ak_array_of_empty_arrays(n):
    content = ak.Array(np.empty(0, dtype=np.float64))
    return ak.unflatten(content, np.zeros(n, dtype=np.int64))


def _next_rc_evt_file(evt_files: list[str | Path], rc_file_state: dict) -> str | Path:
    """Return the next `evt` file, cycling through all files in shuffled order before repeating."""
    if "order" not in rc_file_state:
        order = list(evt_files)
        random.shuffle(order)
        rc_file_state["order"] = order
        rc_file_state["idx"] = 0
        rc_file_state["counts"] = {}

    order = rc_file_state["order"]
    idx = rc_file_state["idx"]

    # if all files were used once, shuffle again and reset index to 0
    if idx >= len(order):
        random.shuffle(order)
        idx = 0

    evt_file = order[idx]
    evt_file_key = str(evt_file)
    reuse_count = rc_file_state["counts"].get(evt_file_key, 0)
    if reuse_count >= 1:
        log.warning(
            "reusing forced-trigger EVT file %s for RC correction (%d previous uses)",
            Path(evt_file).name,
            reuse_count,
        )
    rc_file_state["idx"] = idx + 1
    rc_file_state["counts"][evt_file_key] = reuse_count + 1

    return evt_file


def process_sipm(
    iterator: LH5Iterator,
    optmap_lar: str | Path | OptmapForConvolve,
    sipm: str,
    sipm_uid: int,
    out_file: str | Path,
    runid: str,
    usability: str,
    rc_evt_files: list[str | Path] | None,
    rc_file_state: dict,
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
                map_scaling=MAP_SCALING,
                max_pes_per_hit=MAX_PES_PER_HIT,
            )
        if MAX_PES_PER_HIT > 0:
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
                pe_times_micro, DEFAULT_PHOTOELECTRON_RES
            )

        if TIME_RESOLUTION_NS > 0:
            with perf_block("cluster_photoelectrons()"):
                pe_times, pe_amps = reboost_utils.cluster_photoelectrons(
                    pe_times_micro,
                    pe_amps_micro,
                    TIME_RESOLUTION_NS,
                )
        else:
            pe_times = pe_times_micro
            pe_amps = pe_amps_micro

        # Add random coincidences from forced trigger files
        if rc_evt_files is not None:
            with perf_block("add_random_coincidences()"):
                if usability == "on":
                    # Build a per-chunk RC library, adding more evt files if
                    # one file does not contain enough events for the chunk.
                    rc_parts: list[ak.Array] = []
                    total_rc_events = 0
                    attempts = 0
                    max_attempts = len(rc_evt_files)

                    # Reuse overshoot from previous chunk before loading files.
                    carryover = rc_file_state.get("carryover")
                    if carryover is not None and len(carryover) > 0:
                        carryover_len = len(carryover)
                        n_from_carryover = min(len(chunk), len(carryover))
                        rc_parts.append(carryover[:n_from_carryover])
                        total_rc_events += n_from_carryover
                        if n_from_carryover < len(carryover):
                            rc_file_state["carryover"] = carryover[n_from_carryover:]
                        else:
                            rc_file_state["carryover"] = None
                        remaining_carryover = rc_file_state.get("carryover")
                        remaining_carryover_len = (
                            len(remaining_carryover)
                            if remaining_carryover is not None
                            else 0
                        )
                        log.debug(
                            "used %d residual RC events from carryover; %d remained queued "
                            "(had %d, chunk needs %d total)",
                            n_from_carryover,
                            remaining_carryover_len,
                            carryover_len,
                            len(chunk),
                        )

                    while total_rc_events < len(chunk) and attempts < max_attempts:
                        rc_evt_file = _next_rc_evt_file(rc_evt_files, rc_file_state)
                        n_missing = len(chunk) - total_rc_events
                        part = spms_pars.get_rc_library([rc_evt_file], n_missing)
                        attempts += 1
                        if len(part) == 0:
                            log.warning(
                                "forced-trigger library from %s is empty for this chunk; "
                                "trying next file",
                                Path(rc_evt_file).name,
                            )
                            continue

                        # Take only what is still needed from the beginning of this file.
                        n_take = min(n_missing, len(part))
                        rc_parts.append(part[:n_take])
                        total_rc_events += n_take

                        # Keep the rest for next chunk(s), rather than discarding it.
                        n_left = len(part) - n_take
                        if n_left > 0:
                            rc_file_state["carryover"] = part[n_take:]
                            log.debug(
                                "stored %d residual RC events from %s after taking %d/%d "
                                "for this chunk",
                                n_left,
                                Path(rc_evt_file).name,
                                n_take,
                                len(part),
                            )

                    if total_rc_events == 0:
                        msg = "no random coincidences available from any evt file"
                        raise RuntimeError(msg)

                    chunk_rc_library = (
                        ak.concatenate(rc_parts) if len(rc_parts) > 1 else rc_parts[0]
                    )

                    rc_amps, rc_times = spms_pars.get_sipm_rc_data(
                        chunk_rc_library, sipm, sipm_uid
                    )
                else:
                    rc_amps = _ak_array_of_empty_arrays(len(chunk))
                    rc_times = _ak_array_of_empty_arrays(len(chunk))

        with perf_block("write_chunk()"):
            out_table = reboost_utils.make_output_chunk(lgdo_chunk)

            out_table.add_field(
                "time", VectorOfVectors(pe_times, attrs={"units": "ns"})
            )
            out_table.add_field("energy", VectorOfVectors(pe_amps))
            out_table.add_field("is_saturated", Array(is_saturated))

            if rc_evt_files is not None:
                # TODO: units should be automatically forwarded
                out_table.add_field(
                    "rc_time", VectorOfVectors(rc_times, attrs={"units": "ns"})
                )
                out_table.add_field("rc_energy", VectorOfVectors(rc_amps))

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

    # Lookup forced-trigger evt files for this partition
    if add_random_coincidences:
        msg = "looking up forced trigger files for random coincidences"
        log.debug(msg)
        with perf_block("lookup_rc_files()"):
            evt_tier_name = utils.get_evt_tier_name(l200data)
            rc_evt_files = spms_pars.lookup_evt_files(l200data, runid, evt_tier_name)
            if len(rc_evt_files) == 0:
                msg = "no random coincidences available, cannot continue"
                raise RuntimeError(msg)
    else:
        rc_evt_files = None

    # File-state tracker(s) for iterating through forced-trigger evt files
    # - optmap_per_sipm=True: keep a separate state per SiPM
    # - optmap_per_sipm=False: keep a single shared state ("all")
    rc_file_state = {}

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
                buffer_len=BUFFER_LEN,
            )

        if optmap_per_sipm:
            for sipm in sorted(reboost_utils.get_senstables(geom, "optical")):
                sipm_uid = sens_tables[sipm].uid

                # get the usability
                usability = usabilities[runid][sipm]
                if usability is None:
                    usability = "on"

                msg = f"applying optical map for SiPM {sipm}"
                log.debug(msg)

                process_sipm(
                    _make_iterator(),
                    optmap_lar,
                    sipm,
                    sipm_uid,
                    hit_file,
                    runid,
                    usability,
                    rc_evt_files,
                    rc_file_state.setdefault(sipm, {}),
                )

                print_perf_last()
        else:
            log.debug("applying sum optical map")

            process_sipm(
                _make_iterator(),
                optmap_lar,
                "all",
                geom_meta.uid,
                hit_file,
                runid,
                "on",
                rc_evt_files,
                rc_file_state,
            )


log.debug("building the TCM")
build_tcm(hit_file, hit_file)

with perf_block("move_to_cfs()"):
    move2cfs()

print_perf()
