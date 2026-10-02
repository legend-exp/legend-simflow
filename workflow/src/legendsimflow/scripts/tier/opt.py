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
import logging
from collections.abc import Mapping
from pathlib import Path

import awkward as ak
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
import numpy as np
import pyg4ometry
import pygeomtools
import reboost
import reboost.spms
from dbetto import AttrsDict
from dbetto.utils import load_dict
from lgdo import Array, Table, VectorOfVectors
from lh5 import LH5Iterator
from reboost.optmap.convolve import NumdetStats, OptmapForConvolve
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import metadata as mutils
from legendsimflow import nersc, utils
from legendsimflow import reboost as reboost_utils
from legendsimflow.exceptions import SimflowConfigError
from legendsimflow.metadata import get_tier_settings
from legendsimflow.scripts import log_script_invocation


def resolve_map_scaling(setting: float | Mapping[str, float], sipm: str) -> float:
    """Return the optical-map scaling factor to apply to one SiPM.

    ``optmap_scaling_factor`` is either a scalar applied to every channel, or a
    mapping ``<sipm name> -> <scaling factor>`` for per-channel photon detection
    efficiencies. In the latter case channels missing from the mapping fall back
    to the reserved ``default`` key; without it, a missing channel is an error
    rather than a silent guess.

    Parameters
    ----------
    setting
        The ``optmap_scaling_factor`` value from the ``opt`` tier settings.
    sipm
        Name of the SiPM being processed, or ``"all"`` for the combined map.

    Returns
    -------
    float
    """
    if not isinstance(setting, Mapping):
        return float(setting)

    if sipm in setting:
        return float(setting[sipm])

    if "default" in setting:
        return float(setting["default"])

    msg = (
        f"no optmap_scaling_factor entry for {sipm} and no 'default' key in the "
        "opt tier settings"
    )
    raise KeyError(msg)


@snakemake_compatible(
    mapping={
        "stp_file": "input.stp_file",
        "optmap_lar": lambda snakemake: (
            snakemake.input.optmap_lar[0] if snakemake.input.optmap_lar else None
        ),
        "geom_file": "input.geom",
        "simstat_part_file": "input.simstat_part_file",
        "usability_file": "input.usability",
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
    parser.add_argument(
        "--optmap-lar",
        default=None,
        help="LAr optical map file, needed with light_source: optmap",
    )
    parser.add_argument("--geom-file", required=True, help="GDML geometry file")
    parser.add_argument(
        "--simstat-part-file",
        required=True,
        help="simulation statistics partition file",
    )
    parser.add_argument(
        "--usability-file",
        required=True,
        help="detector usability YAML file",
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
    gdml_file = nersc.dvs_ro(config, args.geom_file)
    log_file = args.log_file
    metadata = config.metadata
    optmap_per_sipm = args.optmap_per_sipm
    scintillator_volume_name = args.scintillator_volume_name
    simstat_part_file = nersc.dvs_ro(config, args.simstat_part_file)
    usability_map = AttrsDict(load_dict(nersc.dvs_ro(config, args.usability_file)))

    opt_file, move2cfs = nersc.make_on_scratch(config, opt_file)

    tier_opt_settings = get_tier_settings(config, "opt")
    optmap_scaling_factor = tier_opt_settings.optmap_scaling_factor
    photoelectron_resolution_sigma = tier_opt_settings.photoelectron_resolution_sigma
    time_resolution_in_ns = tier_opt_settings.time_resolution_in_ns
    buffer_len = tier_opt_settings.buffer_len
    store_expected_pes = tier_opt_settings.get("store_expected_pes", False)
    light_source = tier_opt_settings.get("light_source", "optmap")
    # tracked photons always give one table per SiPM
    max_pes_per_hit = (
        tier_opt_settings.max_pes_per_hit_per_sipm
        if optmap_per_sipm or light_source == "tracked_photons"
        else tier_opt_settings.max_pes_per_hit_combined
    )

    settings_block = f"simprod.config.tier.opt.{config.experiment}.settings"
    if light_source not in ("optmap", "tracked_photons"):
        msg = (
            f"light_source must be 'optmap' or 'tracked_photons', not {light_source!r}"
        )
        raise SimflowConfigError(msg, settings_block)
    if light_source == "optmap":
        if args.optmap_lar is None:
            msg = "--optmap-lar is required with light_source: optmap"
            raise SimflowConfigError(msg, settings_block)
        optmap_lar = nersc.dvs_ro(config, args.optmap_lar)
    else:
        optmap_lar = None

    # setup logging
    log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "tier-opt", parser, args)
    perf_block, print_perf, print_perf_last = reboost.make_profiler()

    # load the geometry and retrieve registered sensitive volume tables
    geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
    sens_tables = pygeomtools.detectors.get_all_senstables(geom)

    # fail early and loudly: without a matching volume the loop below would
    # silently do nothing and the output file would never be created
    scintillators = [
        name
        for name, meta in sens_tables.items()
        if meta.detector_type == "scintillator"
    ]
    if scintillator_volume_name not in scintillators:
        msg = (
            f"scintillator volume {scintillator_volume_name} not found in "
            f"{gdml_file}. scintillator volumes in the geometry: "
            f"{sorted(scintillators)}"
        )
        raise SimflowConfigError(msg, settings_block)

    def process_sipm(
        iterator: LH5Iterator,
        optmap_lar: str | Path | OptmapForConvolve | None,
        sipm: str,
        sipm_uid: int,
        out_file: str | Path,
        runid: str,
        usability: str,
    ) -> None:
        """Write the photoelectrons of one SiPM.

        `iterator` runs over the argon table with ``light_source: optmap`` and over
        the table of the SiPM with ``light_source: tracked_photons``. Each of
        its rows gives one output row.
        """
        if light_source == "optmap":
            with perf_block("load_optmap()"):
                # in per-SiPM mode the map is (re)loaded here, once per call, and
                # released again on return: the per-channel maps are too large to
                # keep all of them resident for the whole job.
                if not isinstance(optmap_lar, OptmapForConvolve):
                    optmap_lar = reboost.spms.load_optmap(optmap_lar, sipm)

            # constant for this SiPM, so resolve it once instead of per chunk
            map_scaling = resolve_map_scaling(optmap_scaling_factor, sipm)
            msg = f"using optical map scaling factor {map_scaling} for {sipm}"
            log.debug(msg)

        total_detected_pe_stats = NumdetStats()

        for lgdo_chunk in iterator:
            chunk = lgdo_chunk.view_as("ak")

            if light_source == "tracked_photons":
                # the SiPM efficiency is applied while tracking: every recorded
                # photon is a photoelectron
                pe_times_micro = ak.sort(chunk.time, axis=-1)
                expected_pes = None
                if max_pes_per_hit > 0:
                    # as with the map: saturated once the cap is reached
                    is_saturated = ak.to_numpy(
                        ak.num(pe_times_micro) >= max_pes_per_hit
                    )
                    pe_times_micro = pe_times_micro[:, :max_pes_per_hit]
                else:
                    is_saturated = np.full(len(chunk), fill_value=False, dtype=np.bool_)
            else:
                with perf_block("emitted_scintillation_photons()"):
                    scint_ph = reboost.spms.emitted_scintillation_photons(
                        chunk.edep, chunk.particle, "lar"
                    )

                with perf_block("number_of_detected_photoelectrons()"):
                    # return_stats=True also silences the per-chunk warnings of
                    # reboost about steps outside the map
                    *_output, _detected_pe_stats = (
                        reboost.spms.number_of_detected_photoelectrons(
                            chunk.xloc,
                            chunk.yloc,
                            chunk.zloc,
                            scint_ph,
                            optmap_lar,
                            sipm,
                            map_scaling=map_scaling,
                            max_pes_per_hit=max_pes_per_hit,
                            return_pes_expectation_value=store_expected_pes,
                            return_stats=True,
                        )
                    )
                total_detected_pe_stats += _detected_pe_stats

                # reboost appends the expectation after the other outputs
                if store_expected_pes:
                    *_output, expected_pes = _output
                else:
                    expected_pes = None
                if max_pes_per_hit > 0:
                    nr_pe, is_saturated = _output
                else:
                    (nr_pe,) = _output
                    is_saturated = np.full(len(chunk), fill_value=False, dtype=np.bool_)

                with perf_block("photoelectron_times()"):
                    pe_times_micro = reboost.spms.photoelectron_times(
                        nr_pe, chunk.particle, chunk.time, "lar"
                    )

                    # the photoelectron_times() processor does not guarantee time
                    # ordering
                    pe_times_micro = ak.sort(pe_times_micro, axis=-1)

            with perf_block("photoelectron_resolution()"):
                pe_amps_micro = reboost.spms.smear_photoelectrons(
                    pe_times_micro, photoelectron_resolution_sigma
                )

            if time_resolution_in_ns > 0:
                with perf_block("cluster_photoelectrons()"):
                    pe_times, pe_amps = reboost.spms.cluster_photoelectrons(
                        pe_times_micro,
                        pe_amps_micro,
                        time_resolution_in_ns,
                    )
            else:
                pe_times = pe_times_micro
                pe_amps = pe_amps_micro

            with perf_block("write_hit_table_chunk()"):
                out_table = reboost.init_hit_table(lgdo_chunk)

                # relative to the hit t0, subtracted in float64
                pe_times = pe_times - out_table.t0.view_as("ak")
                out_table.add_field(
                    "dt",
                    VectorOfVectors(
                        ak.values_astype(pe_times, np.float32), attrs={"units": "ns"}
                    ),
                )
                out_table.add_field(
                    "energy", VectorOfVectors(ak.values_astype(pe_amps, np.float32))
                )
                out_table.add_field("is_saturated", Array(is_saturated))
                if expected_pes is not None:
                    out_table.add_field(
                        "expected_pes",
                        Array(np.asarray(expected_pes, dtype=np.float32)),
                    )

                _, period, run, _ = mutils.parse_runid(runid)
                field_vals = [period, run, mutils.encode_usability(usability)]
                for i, field in enumerate(["period", "run", "usability"]):
                    out_table.add_field(
                        field,
                        Array(np.full(shape=len(chunk), fill_value=field_vals[i])),
                    )

                reboost.write_hit_table_chunk(
                    out_table,
                    "hit/" + ("spms" if sipm == "all" else sipm),
                    out_file,
                    uid=sipm_uid,
                )

        if light_source == "tracked_photons":
            return

        tot = total_detected_pe_stats.energy_looped
        for counts, where in (
            (
                total_detected_pe_stats.energy_oob
                + total_detected_pe_stats.energy_no_stats,
                "outside the map",
            ),
            (
                total_detected_pe_stats.energy_zero,
                "in bins with zero detection probability",
            ),
        ):
            pct = 100 * counts / tot if tot > 0 else 0.0
            msg = (
                f"optical map {sipm}: {pct:.2f}% of the energy deposited "
                f"in {scintillator_volume_name} is {where}"
            )
            log.log(logging.WARNING if pct > 1 else logging.DEBUG, msg)

    partitions = load_dict(simstat_part_file)[f"job_{jobid}"]

    # load TCM, to be used to chunk the event statistics according to the run partitioning
    msg = "loading TCM"
    log.debug(msg)
    tcm = lh5.read_as("tcm", stp_file, library="ak")

    # in combined mode there is a single map, so pre-load it once for a little
    # speed up. in per-SiPM mode the maps are loaded on demand in process_sipm()
    if light_source == "optmap" and not optmap_per_sipm:
        optmap_lar = reboost.spms.load_optmap(optmap_lar, "all")

    sipms = sorted(reboost_utils.get_senstables(geom, "optical"))

    # loop over the partitions for this file
    for runid_idx, (runid, evt_idx_range) in enumerate(partitions.items()):
        msg = (
            f"processing partition corresponding to {runid} "
            f"[{runid_idx + 1}/{len(partitions)}], event range {evt_idx_range}"
        )
        log.info(msg)

        if light_source == "tracked_photons":
            # as in the hit tier, each row of a SiPM table gives one output row
            for sipm in sipms:
                sipm_uid = sens_tables[sipm].uid
                i_start, n_entries = reboost.get_rows_in_event_range(
                    tcm, sipm_uid, *evt_idx_range
                )
                if n_entries == 0:
                    continue

                usability = usability_map[runid].get(sipm)
                if usability is None:
                    msg = f"usability not found for {sipm} in {runid}, defaulting to on"
                    log.warning(msg)
                    usability = "on"

                msg = f"processing the tracked photons of SiPM {sipm}"
                log.debug(msg)

                process_sipm(
                    LH5Iterator(
                        stp_file,
                        f"stp/{sipm}",
                        i_start=i_start,
                        n_entries=n_entries,
                        buffer_len=buffer_len,
                        field_mask=["evtid", "t0", "time"],
                    ),
                    None,
                    sipm,
                    sipm_uid,
                    opt_file,
                    runid,
                    usability,
                )
                print_perf_last()
        else:
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

                if f"stp/{det_name}" not in lh5.ls(stp_file, "stp/*"):
                    msg = (
                        f"detector {det_name} not found in {stp_file}. "
                        "possibly because it was not read-out or there were no hits recorded"
                    )
                    log.warning(msg)
                    continue

                log.info("processing the 'lar' scintillator table...")

                msg = "looking for indices of hit table rows to read..."
                log.debug(msg)
                i_start, n_entries = reboost.get_rows_in_event_range(
                    tcm, geom_meta.uid, *evt_idx_range
                )

                def _make_iterator(
                    det_name=det_name, i_start=i_start, n_entries=n_entries
                ):
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
                        usability = usability_map[runid].get(sipm)
                        if usability is None:
                            msg = f"usability not found for {sipm} in {runid}, defaulting to on"
                            log.warning(msg)
                            usability = "on"

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

    # a SiPM without photons still gets a table: the evt tier takes the list of
    # channels from the opt file, which must not change from job to job
    if light_source == "tracked_photons":
        no_photons = Table(
            {
                "evtid": Array(np.empty(0, dtype=np.int64)),
                "t0": Array(np.empty(0), attrs={"units": "ns"}),
                "time": VectorOfVectors(
                    flattened_data=Array(np.empty(0)),
                    cumulative_length=Array(np.empty(0, dtype=np.int64)),
                    attrs={"units": "ns"},
                ),
            }
        )
        tables = lh5.ls(opt_file, "hit/*") if Path(opt_file).exists() else []
        for sipm in sipms:
            if f"hit/{sipm}" in tables:
                continue
            process_sipm(
                [no_photons], None, sipm, sens_tables[sipm].uid, opt_file, runid, "on"
            )

    log.debug("building the TCM")
    reboost.build_remage_tcm(opt_file, opt_file)

    with perf_block("move_to_cfs()"):
        move2cfs()

    print_perf()


if __name__ == "__main__":
    main()
