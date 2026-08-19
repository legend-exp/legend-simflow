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

import awkward as ak
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
import numpy as np
import reboost.pmts as reboost_pmts
from dbetto import AttrsDict
from dbetto.utils import load_dict
from lgdo import Array, Scalar, Struct, Table, VectorOfVectors
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import nersc, spms_pars, utils
from legendsimflow import reboost as reboost_utils
from legendsimflow.awkward import ak_isin
from legendsimflow.metadata import (
    encode_psd_usability,
    encode_usability,
    get_tier_settings,
    parse_runid,
)
from legendsimflow.profile import make_profiler
from legendsimflow.scripts import log_script_invocation
from legendsimflow.tcm import merge_stp_n_opt_tcms_to_lh5

OFF = encode_usability("off")
ON = encode_usability("on")
VALID_PSD = encode_psd_usability("valid")


@snakemake_compatible(
    mapping={
        "stp_file": "input.stp_file",
        "opt_file": lambda snakemake: (
            snakemake.input.opt_file[0] if snakemake.input.opt_file else None
        ),
        "hit_file": lambda snakemake: (
            snakemake.input.hit_file[0] if snakemake.input.hit_file else None
        ),
        "simstat_part_file": "input.simstat_part_file",
        "usability_file": "input.usability",
        "jobid": "wildcards.jobid",
        "evt_file": "output[0]",
        "log_file": "log[0]",
        "add_random_coincidences": "params.add_random_coincidences",
        "skip_opt": "params.skip_opt",
        "skip_hit": "params.skip_hit",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(description="Build the evt tier.")
    parser.add_argument("--stp-file", required=True, help="input stp tier file")
    parser.add_argument(
        "--opt-file", required=False, default=None, help="input opt tier file"
    )
    parser.add_argument(
        "--hit-file", required=False, default=None, help="input hit tier file"
    )
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
    parser.add_argument("--evt-file", required=True, help="output evt tier file")
    parser.add_argument("--log-file", default=None, help="log file")
    parser.add_argument(
        "--add-random-coincidences",
        action="store_true",
        default=False,
        help="add random-coincidence SiPM data from real evt files",
    )
    parser.add_argument(
        "--skip-opt",
        action="store_true",
        default=False,
        help="skip the opt (SiPM) tier; no spms table will be produced",
    )
    parser.add_argument(
        "--skip-hit",
        action="store_true",
        default=False,
        help="skip the hit (HPGe) tier; no geds table will be produced",
    )
    parser.add_argument(
        "--simflow-config",
        "--config",
        dest="simflow_config",
        required=True,
        help="simflow config YAML path",
    )
    args = parser.parse_args()

    skip_opt = args.skip_opt
    skip_hit = args.skip_hit

    if skip_opt and skip_hit:
        parser.error("cannot set both --skip-opt and --skip-hit")
    if not skip_opt and args.opt_file is None:
        parser.error("--opt-file is required unless --skip-opt is set")
    if not skip_hit and args.hit_file is None:
        parser.error("--hit-file is required unless --skip-hit is set")

    config = utils.init_simflow_context(args.simflow_config, workflow=None).config

    stp_file = nersc.dvs_ro(config, args.stp_file)
    hit_file = {}
    if not skip_opt:
        hit_file["opt"] = nersc.dvs_ro(config, args.opt_file)
    if not skip_hit:
        hit_file["hit"] = nersc.dvs_ro(config, args.hit_file)
    evt_file = args.evt_file
    log_file = args.log_file
    metadata = config.metadata
    tier_evt_settings = get_tier_settings(config, "evt")
    geds_energy_thr_kev = tier_evt_settings.geds_energy_thr_kev
    spms_energy_thr_pe = tier_evt_settings.spms_energy_thr_pe
    skip_muonveto = tier_evt_settings.skip_muonveto
    muonveto_trigger_groups = tier_evt_settings.muonveto_trigger_groups
    muonveto_trigger_timegate = tier_evt_settings.muonveto_trigger_timegate
    muonveto_trigger_deadtime = tier_evt_settings.muonveto_trigger_deadtime
    muonveto_trace_length = tier_evt_settings.muonveto_trace_length
    muonveto_time_per_sample = tier_evt_settings.muonveto_time_per_sample
    muonveto_trigger_position = tier_evt_settings.muonveto_trigger_position
    muonveto_coincidence_window = tier_evt_settings.muonveto_coincidence_window
    buffer_len = tier_evt_settings.buffer_len
    simstat_part_file = nersc.dvs_ro(config, args.simstat_part_file)
    add_random_coincidences = args.add_random_coincidences
    l200data = config.paths.get("l200data", None)
    usability_map = AttrsDict(load_dict(nersc.dvs_ro(config, args.usability_file)))

    # get the psd settings
    tier_hit_settings = get_tier_settings(config, "hit")
    simulate_psd = tier_hit_settings.get("simulate_psd", True)
    simulate_psd_with_psl = tier_hit_settings.get("simulate_psd_with_psl", False)

    evt_file, move2cfs = nersc.make_on_scratch(config, evt_file)

    # setup logging
    log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "tier-evt", parser, args)
    perf_block, print_stats, print_stats_since_last = make_profiler()

    log.info("merging hit and opt TCMs")
    with perf_block("merge_tcms()"):
        if skip_opt:
            # No SiPM data: stream STP TCM directly into evt file
            wo = "write_safe"
            for chunk in lh5.LH5Iterator(str(stp_file), "tcm", buffer_len=buffer_len):
                lh5.write(chunk, "tcm", str(evt_file), wo_mode=wo)
                wo = "append"
        else:
            scintillator_uid = next(
                uid
                for uid, name in reboost_utils.get_remage_detector_uids(
                    stp_file
                ).items()
                if name == "liquid_argon"
            )

            merge_stp_n_opt_tcms_to_lh5(
                stp_file,
                hit_file["opt"],
                evt_file,
                scintillator_uid=scintillator_uid,
                buffer_len=buffer_len,
            )

    # test that the evt tcm has the same amount of rows as the stp tcm
    if lh5.read_n_rows("tcm", stp_file) != lh5.read_n_rows("tcm", evt_file):
        msg = (
            "stp and evt tcm should have same number of rows not "
            f"stp={lh5.read_n_rows('tcm', stp_file)}, "
            f"evt={lh5.read_n_rows('tcm', evt_file)}"
        )
        if not skip_hit:
            msg += f", hit={lh5.read_n_rows('tcm', hit_file['hit'])}"
        if not skip_opt:
            msg += f", opt={lh5.read_n_rows('tcm', hit_file['opt'])}"
        raise ValueError(msg)

    # get the mapping of detector name to uid
    # NOTE: we check on disk because we are not sure which tables were processed in
    # the hit tiers
    det2uid = {}
    for tier in ("opt", "hit", "stp"):
        if (tier == "opt" and skip_opt) or (tier == "hit" and skip_hit):
            det2uid[tier] = {}
            continue
        det2uid[tier] = {
            name: uid
            for uid, name in reboost_utils.get_remage_detector_uids(
                hit_file[tier], lh5_table="hit"
            ).items()
            if (tier != "stp") or ("pmt" in name)  # for the stp tier only map the pmts
        }
        msg = f"found mapping name -> uid ({tier} tier): {det2uid[tier]}"
        log.debug(msg)

    # little helper to simplify the code below
    # TODO: move/fix in reboost
    def _read_hits(tcm_ak, tier, field):
        if not det2uid[tier]:
            return ak.Array([[] for _ in range(len(ak.num(tcm_ak[tier].row_in_table)))])

        msg = f"loading {field=} data from {tier=} (file {hit_file[tier]})"
        log.debug(msg)

        tcm = tcm_ak[tier]
        tcm_flat = ak.Array({k: ak.flatten(tcm[k]) for k in tcm.fields})

        data_flat = []
        tcm_rows = []

        # for un-flattening at the end
        counts = ak.num(tcm.row_in_table)

        for tab_name, key in det2uid[tier].items():
            mask = tcm_flat.table_key == key

            with perf_block("_read_hits()_tcm_filter"):
                rows = np.sort(tcm_flat.row_in_table[mask].to_numpy())
                tcm_rows.append(np.where(mask)[0].to_numpy())

            with perf_block("_read_hits()_lh5.read()"):
                # check if we can just use the start_row / n_rows arguments
                # to read. this seems to be faster than using the idx argument
                # TODO: check/fix in legend-lh5io
                if len(rows) >= 2 and np.all(rows == np.arange(rows[0], rows[-1] + 1)):
                    data_ch = lh5.read(
                        f"hit/{tab_name}/{field}",
                        hit_file[tier],
                        start_row=rows[0],
                        n_rows=len(rows),
                    )
                else:
                    msg = (
                        "unexpected: hit rows indices are != range(rows[0], rows[-1]+1). "
                        "falling back to using lh5.read(..., idx=rows)"
                    )
                    log.warning(msg)
                    data_ch = lh5.read(
                        f"hit/{tab_name}/{field}", hit_file[tier], idx=rows
                    )

            units = data_ch.attrs.get("units", None)
            data_ch = data_ch.view_as("ak")

            data_flat.append(data_ch)

        tcm_rows_concat = np.concatenate(tcm_rows)
        data_flat_concat = ak.concatenate(data_flat)[np.argsort(tcm_rows_concat)]

        data_unflat = ak.unflatten(data_flat_concat, counts)

        if units is not None:
            return ak.with_parameter(data_unflat, "units", units)
        return data_unflat

    partitions = load_dict(simstat_part_file)[f"job_{args.jobid}"]

    # use write_safe on the first chunk to catch stale data from a failed retry
    evt_wo_mode = "write_safe"

    log.info("begin iterating over TCM")
    for runid_idx, (runid, evt_idx_range) in enumerate(partitions.items()):
        msg = (
            f"processing partition corresponding to {runid} "
            f"[{runid_idx + 1}/{len(partitions)}], event range {evt_idx_range}"
        )
        log.info(msg)

        evt_start, evt_end = evt_idx_range
        # evt_idx_range is [start, end] inclusive
        n_entries = evt_end - evt_start + 1

        # canonical non-OFF SiPM channel UIDs for this run, in ascending order.
        # used to pad events with no edep in LAr photons
        on_spms_uids = (
            sorted(
                uid
                for det_name, uid in det2uid["opt"].items()
                if usability_map[runid].get(det_name, "on") != "off"
            )
            if not skip_opt
            else []
        )

        if add_random_coincidences:
            msg = "looking up forced trigger files for random coincidences"
            log.debug(msg)
            with perf_block("lookup_l200data_evts_for_rc()"):
                evt_tier_name = utils.get_evt_tier_name(l200data)
                rc_evt_files = sorted(
                    spms_pars.lookup_evt_files(l200data, runid, evt_tier_name)
                )
                if not rc_evt_files:
                    msg = "no RC evt files found for random coincidences"
                    raise RuntimeError(msg)

                rc_index_lookup = spms_pars.build_rc_evt_index_lookup(rc_evt_files)
            # state is reset per partition so RC events are drawn independently
            # for each run slice
            rc_file_state: dict = {}

        # iterate over the unified tcm for this partition; an empty partition
        # (n_entries=0) produces no chunks and is silently skipped
        # NOTE: open mode is append because we will write to the same file
        it = lh5.LH5Iterator(
            str(evt_file),
            "tcm",
            i_start=evt_start,
            n_entries=n_entries,
            buffer_len=buffer_len,
            h5py_open_mode="a",
        )
        for chunk in it:
            unified_tcm = chunk.view_as("ak")
            out_table = Table(size=len(unified_tcm))

            # split the unified TCM in two, one for each tier. in this way we will be
            # able to read data from each tier
            tcm = {}
            for tier in ("opt", "hit", "stp"):
                mask = ak_isin(unified_tcm.table_key, det2uid[tier].values())
                tcm[tier] = unified_tcm[mask]

            # trigger table
            # -------------
            out_table.add_field("trigger", Table(size=len(unified_tcm)))

            # global fields that are constant over the full events
            # take them from the hit tier if available, otherwise from opt
            trigger_source = "opt" if skip_hit else "hit"
            for constant_field in ["run", "period", "evtid"]:
                if skip_hit and constant_field in ("run", "period"):
                    # opt TCM can be empty for events with no LAr deposit — reading
                    # from opt would give NaN. Source run/period from the partition
                    # runid instead, where they are always known.
                    _, _period, _run, _ = parse_runid(runid)
                    val = _period if constant_field == "period" else _run
                    out_table.add_field(
                        f"trigger/{constant_field}",
                        Array(np.full(len(unified_tcm), val, dtype=np.float64)),
                    )
                    continue

                data = _read_hits(tcm, trigger_source, constant_field)

                # sanity check
                assert len(data) == len(tcm[trigger_source])

                # replace the awkward missing values with NaN for LH5 compatibility
                data = ak.fill_none(ak.firsts(data, axis=-1), np.nan)
                out_table.add_field(f"trigger/{constant_field}", Array(data))

                # Store the event-ids for comparison with muons later.
                if constant_field == "evtid":
                    event_ids = data.to_numpy()

            timestamp = _read_hits(tcm, trigger_source, "t0")
            timestamp = ak.fill_none(ak.firsts(timestamp, axis=-1), np.nan)
            out_table.add_field(
                "trigger/timestamp",
                Array(np.asarray(timestamp, dtype=np.float32), attrs={"units": "ns"}),
            )
            # HPGe table
            # ----------
            if not skip_hit:
                out_table.add_field("geds", Table(size=len(unified_tcm)))

                # first read usability and energy
                usability = _read_hits(tcm, "hit", "usability")
                psd_usability = _read_hits(tcm, "hit", "psd_usability")
                is_valid_sim = _read_hits(tcm, "hit", "is_valid_sim")
                energy = _read_hits(tcm, "hit", "energy")

                # we want to only store hits from events in ON and AC detectors and above
                # our energy threshold
                hitsel = (usability != OFF) & (energy > geds_energy_thr_kev)

                # we want to still be able to know which detectors are ON (and not AC)
                out_table.add_field(
                    "geds/is_good_channel", VectorOfVectors(usability[hitsel] == ON)
                )
                out_table.add_field(
                    "geds/energy",
                    VectorOfVectors(
                        ak.values_astype(energy[hitsel], np.float32),
                        attrs={"units": "keV"},
                    ),
                )
                # NOTE: the energy sum does not include AC detectors
                out_table.add_field(
                    "geds/energy_sum",
                    Array(
                        np.asarray(
                            ak.sum(energy[hitsel & (usability == ON)], axis=-1),
                            dtype=np.float32,
                        ),
                        attrs={"units": "keV"},
                    ),
                )

                # fields to identify detectors and lookup stuff in the lower tiers
                out_table.add_field(
                    "geds/rawid", VectorOfVectors(tcm["hit"].table_key[hitsel])
                )
                out_table.add_field(
                    "geds/hit_idx", VectorOfVectors(tcm["hit"].row_in_table[hitsel])
                )

                # PSD subtable: common flags live directly under psd, while
                # method-specific fields go into dedicated subtables
                if simulate_psd or simulate_psd_with_psl:
                    out_table.add_field("geds/psd", Table(size=len(unified_tcm)))
                    out_table.add_field(
                        "geds/psd/is_good",
                        VectorOfVectors(psd_usability[hitsel] == VALID_PSD),
                    )
                    # whether the inputs needed to model PSD for this hit are
                    # valid (computed in the hit tier)
                    out_table.add_field(
                        "geds/psd/is_valid_sim",
                        VectorOfVectors(is_valid_sim[hitsel]),
                    )

                # single-template A/E PSD
                if simulate_psd:
                    out_table.add_field(
                        "geds/psd/single_temp", Table(size=len(unified_tcm))
                    )

                    aoe = _read_hits(tcm, "hit", "psd/single_temp/aoe")
                    out_table.add_field(
                        "geds/psd/single_temp/aoe",
                        VectorOfVectors(ak.values_astype(aoe[hitsel], np.float32)),
                    )
                    drift_time = _read_hits(
                        tcm, "hit", "psd/single_temp/drift_time_amax"
                    )
                    out_table.add_field(
                        "geds/psd/single_temp/drift_time_amax",
                        VectorOfVectors(
                            ak.values_astype(drift_time[hitsel], np.float32)
                        ),
                    )
                    aoe_corr = _read_hits(tcm, "hit", "psd/single_temp/aoe_corr")
                    out_table.add_field(
                        "geds/psd/single_temp/aoe_corr",
                        VectorOfVectors(ak.values_astype(aoe_corr[hitsel], np.float32)),
                    )
                    out_table.add_field(
                        "geds/psd/single_temp/has_aoe",
                        VectorOfVectors(~np.isnan(aoe[hitsel])),
                    )

                    is_ss = _read_hits(tcm, "hit", "psd/single_temp/is_single_site")
                    out_table.add_field(
                        "geds/psd/single_temp/is_single_site",
                        VectorOfVectors(is_ss[hitsel]),
                    )

                # pulse-shape-library (PSL) based PSD
                if simulate_psd_with_psl:
                    out_table.add_field(
                        "geds/psd/pulse_lib", Table(size=len(unified_tcm))
                    )

                    aoe_corr = _read_hits(tcm, "hit", "psd/pulse_lib/aoe_corr")
                    out_table.add_field(
                        "geds/psd/pulse_lib/aoe_corr",
                        VectorOfVectors(ak.values_astype(aoe_corr[hitsel], np.float32)),
                    )

                    aoe = _read_hits(tcm, "hit", "psd/pulse_lib/aoe")

                    out_table.add_field(
                        "geds/psd/pulse_lib/aoe",
                        VectorOfVectors(ak.values_astype(aoe[hitsel], np.float32)),
                    )
                    out_table.add_field(
                        "geds/psd/pulse_lib/has_aoe",
                        VectorOfVectors(~np.isnan(aoe[hitsel])),
                    )
                    drift_time = _read_hits(tcm, "hit", "psd/pulse_lib/drift_time_amax")
                    out_table.add_field(
                        "geds/psd/pulse_lib/drift_time_amax",
                        VectorOfVectors(
                            ak.values_astype(drift_time[hitsel], np.float32)
                        ),
                    )
                    is_ss = _read_hits(tcm, "hit", "psd/pulse_lib/is_single_site")
                    out_table.add_field(
                        "geds/psd/pulse_lib/is_single_site",
                        VectorOfVectors(is_ss[hitsel]),
                    )
                    is_bb_like = _read_hits(tcm, "hit", "psd/pulse_lib/is_bb_like")
                    out_table.add_field(
                        "geds/psd/pulse_lib/is_bb_like",
                        VectorOfVectors(is_bb_like[hitsel]),
                    )
                    is_high_aoe = _read_hits(tcm, "hit", "psd/pulse_lib/is_high_aoe")
                    out_table.add_field(
                        "geds/psd/pulse_lib/is_high_aoe",
                        VectorOfVectors(is_high_aoe[hitsel]),
                    )
                # compute multiplicity
                geds_multiplicity = ak.sum(hitsel, axis=-1)
                out_table.add_field("geds/multiplicity", Array(geds_multiplicity))
            else:
                geds_multiplicity = ak.Array(np.zeros(len(unified_tcm), dtype=np.int64))

            # SiPM table
            # ----------
            if not skip_opt:
                out_table.add_field("spms", Table(size=len(unified_tcm)))

                # also here, we exclude the non usable channels. this is in line with what
                # done in the evt tier in pygama
                usability = _read_hits(tcm, "opt", "usability")
                energy = _read_hits(tcm, "opt", "energy")
                chansel = usability != OFF
                # we also discard all pulses with amplitude below threshold
                pesel = energy > spms_energy_thr_pe

                # in simulation the opt TCM does not record events for which there is
                # no energy in LAr. This means that in the unified TCM these events
                # will be characterized by spms empty arrays.  pad those events with
                # the canonical non-OFF channel list (empty PE arrays) to match the
                # real-data convention where all non-OFF channels are always present.
                # NOTE: the on_spms_uids ordering must match the ordering
                # used by non-empty events (i.e. the TCM ordering). Currently both
                # are ascending by UID.
                n_events = len(unified_tcm)
                is_empty_opt = ak.num(tcm["opt"].table_key) == 0
                rawid = ak.Array([on_spms_uids] * n_events)

                # rawid is the same canonical list for every event (non-empty events
                # already carry all non-OFF channels in ascending UID order)
                out_table.add_field("spms/rawid", VectorOfVectors(rawid))

                energy_sel = energy[pesel][chansel]
                # fill in empty arrays for events with no LAr edep
                empty_energy = ak.Array([[[] for _ in on_spms_uids]] * n_events)
                energy_sel = ak.where(is_empty_opt, empty_energy, energy_sel)
                out_table.add_field(
                    "spms/energy",
                    VectorOfVectors(ak.values_astype(energy_sel, np.float32)),
                )

                is_saturated = _read_hits(tcm, "opt", "is_saturated")
                is_saturated_sel = is_saturated[chansel]
                # fill in Falses for events with no LAr edep
                empty_is_saturated = ak.Array(
                    [[False for _ in on_spms_uids]] * n_events
                )
                is_saturated_sel = ak.where(
                    is_empty_opt, empty_is_saturated, is_saturated_sel
                )
                out_table.add_field(
                    "spms/is_saturated", VectorOfVectors(is_saturated_sel)
                )

                hit_idx = tcm["opt"].row_in_table[chansel]
                # fill in -1 hit index for events with no LAr edep
                empty_hit_idx = ak.Array([[-1 for _ in on_spms_uids]] * n_events)
                hit_idx = ak.where(is_empty_opt, empty_hit_idx, hit_idx)
                out_table.add_field("spms/hit_idx", VectorOfVectors(hit_idx))

                time = _read_hits(tcm, "opt", "time")
                time_sel = time[pesel][chansel]
                # fill in empty arrays for events with no LAr edep
                empty_time = ak.Array([[[] for _ in on_spms_uids]] * n_events)
                time_sel = ak.where(is_empty_opt, empty_time, time_sel)
                out_table.add_field(
                    "spms/time",
                    VectorOfVectors(
                        ak.values_astype(time_sel, np.float32), attrs={"units": "ns"}
                    ),
                )

                if add_random_coincidences:
                    with perf_block("get_chunk_rc_data()"):
                        rc_chunk = spms_pars.get_chunk_rc_data(
                            [str(f) for f in rc_evt_files],
                            rc_file_state,
                            len(unified_tcm),
                            rc_index_lookup,
                        )
                    # FIXME: this assertion fails because we haven't thought about
                    # cases when there is a DAQ recabling without hardware changes.
                    # right now this fails with p18 because SiPMs were recabled.
                    #
                    # assert rawid alignment: RC and simulation must use the same
                    # channel ordering (both are ascending by UID)
                    # assert ak.to_list(rc_chunk.rawid[0]) == on_spms_uids, (
                    #     "RC rawid does not match simulation spms/rawid: "
                    #     f"{rc_chunk.rawid[0].to_list()} != {on_spms_uids}"
                    # )
                    out_table.add_field(
                        "spms/rc_energy",
                        VectorOfVectors(ak.values_astype(rc_chunk.npe, np.float32)),
                    )
                    out_table.add_field(
                        "spms/rc_time",
                        VectorOfVectors(
                            ak.values_astype(rc_chunk.t0, np.float32),
                            attrs={"units": "ns"},
                        ),
                    )

                # total amount of light per event
                energy_sum = ak.sum(ak.sum(energy[pesel][chansel], axis=-1), axis=-1)
                out_table.add_field(
                    "spms/energy_sum",
                    Array(np.asarray(energy_sum, dtype=np.float32)),
                )

                # how many channels saw some light
                spms_multiplicity = ak.sum(ak.any(chansel & pesel, axis=-1), axis=-1)
                out_table.add_field("spms/multiplicity", Array(spms_multiplicity))
            else:
                energy_sum = ak.Array(np.zeros(len(unified_tcm), dtype=np.float32))
                spms_multiplicity = ak.Array(np.zeros(len(unified_tcm), dtype=np.int64))

            # muonveto table
            # ------------------
            if not skip_muonveto:
                out_table.add_field("muon", Table(size=len(unified_tcm)))
                # TODO : Figure out a way to get the usable pmts of this run from somewhere.
                # I do not think the metadata has this information. Especially as PMTs sometimes break in the middle of a run.
                usable_pmts = list(det2uid[tier].values())

                muonveto_times = _read_hits(tcm, "stp", "time")
                # Group the times with only usable detectors and pad empty events to align the event axis.
                grouped = reboost_pmts.group_by_detector(
                    muonveto_times, tcm["stp"].table_key, usable_pmts
                )
                # TODO : Derive the trigger groups from the usable PMTs and their position (which should be in the metadata)
                triggers = reboost_pmts.build_hardware_triggers(
                    grouped,
                    trigger_groups=muonveto_trigger_groups,
                    timegate=muonveto_trigger_timegate,
                    trigger_deadtime=muonveto_trigger_deadtime,
                )[0]
                hits = reboost_pmts.build_hits(
                    grouped,
                    triggers,
                    trace_length=muonveto_trace_length,
                    time_per_sample=muonveto_time_per_sample,
                    trigger_position=muonveto_trigger_position,
                )
                del (
                    grouped
                )  # grouped consumes a lot of memory and is not needed anymore.
                integral_light = ak.sum(hits, axis=0)
                multiplicity = ak.sum(hits > 0, axis=0)

                # Handle a rare edge case:
                # The tcm can split the same eventid into multiple rows. Make sure every row has all triggers for one G4Event.
                # If this turns out to be not so rare, this part might need to be optimized.
                if len(np.unique(event_ids)) != len(event_ids):
                    unique_ids, counts = np.unique(event_ids, return_counts=True)
                    duplicated_ids = set(unique_ids[counts > 1])
                    triggers_list = ak.to_list(triggers)
                    integral_light_list = ak.to_list(integral_light)
                    multiplicity_list = ak.to_list(multiplicity)
                    for eid in duplicated_ids:
                        idxs = np.where(event_ids == eid)[0]
                        all_triggers = ak.to_list(ak.flatten(triggers[idxs]))
                        all_integral_light = ak.to_list(
                            ak.flatten(integral_light[idxs])
                        )
                        all_multiplicity = ak.to_list(ak.flatten(multiplicity[idxs]))

                        for i in idxs:
                            triggers_list[i] = all_triggers
                            integral_light_list[i] = all_integral_light
                            multiplicity_list[i] = all_multiplicity

                    triggers = ak.Array(triggers_list)
                    integral_light = ak.Array(integral_light_list)
                    multiplicity = ak.Array(multiplicity_list)

                # This is the exact condition from the legend-dataflow-config, except the noisy cut.
                # Might be reasonable to move this to the config? But the LAr cut below is also hard-coded so i did it here...
                software_cut = (
                    (
                        (multiplicity >= 11) | (integral_light >= 7)
                    )  # low-multiplicity cut
                    & ((multiplicity >= 13) | (integral_light <= 98))  # flasher cut
                )
                # software_cut still is a jagged awkward array, as there might be multiple muonveto triggers per tcm event.

                # sanity check
                assert len(integral_light) == len(event_ids)
                assert len(multiplicity) == len(event_ids)

                dt = triggers - timestamp[:, None]
                # muonveto trigger happened after/before the event trigger, but within the coincidence window
                coincidence_mask = (dt >= muonveto_coincidence_window[0]) & (
                    dt <= muonveto_coincidence_window[1]
                )

                # veto only events that are in coincidence with a muonveto trigger and pass the software cut
                veto_mask = coincidence_mask & software_cut
                # If any of the muonveto triggers is in coincidence with the event trigger, we veto the event.
                muonveto_offline = np.any(veto_mask, axis=1)
                muon_flag = np.any(coincidence_mask, axis=1)

                out_table.add_field(
                    "muon/timestamp",
                    VectorOfVectors(
                        ak.values_astype(triggers[coincidence_mask], np.float32),
                        attrs={"units": "ns"},
                    ),
                )
                out_table.add_field(
                    "muon/energy_sum",
                    VectorOfVectors(
                        ak.values_astype(integral_light[coincidence_mask], np.float32)
                    ),
                )
                out_table.add_field(
                    "muon/multiplicity",
                    VectorOfVectors(
                        ak.values_astype(multiplicity[coincidence_mask], np.int32)
                    ),
                )

            # coincidences table
            # ------------------
            out_table.add_field("coincident", Table(size=len(unified_tcm)))

            # is there a signal in the HPGe array?
            if not skip_hit:
                out_table.add_field("coincident/geds", Array(geds_multiplicity > 0))

            # is there a signal in the LAr instrumentation?
            if not skip_opt:
                lar_veto = (spms_multiplicity >= 4) | (energy_sum >= 4)
                out_table.add_field("coincident/spms", Array(lar_veto))

            if not skip_muonveto:
                # This is the "muon" field in the data, which is the aux channel. It is True if there was a muonveto hardware trigger in coincidence
                out_table.add_field("coincident/muon", Array(muon_flag))
                # This represents the "muon_offline" field in the data. True if the coincident muontrigger passed the software cut.
                out_table.add_field("coincident/muon_offline", Array(muonveto_offline))

            # now write down
            with perf_block("write_chunk()"):
                lh5.write(out_table, "evt", evt_file, wo_mode=evt_wo_mode)
                evt_wo_mode = "append"

        print_stats_since_last()

    # always written (empty when both skip_hit and skip_opt) so cvt can rely on
    # a stable schema. Reboost UIDs are disjoint between systems so the union
    # is unambiguous.
    detector_uids = Struct(
        {
            name: Scalar(int(uid))
            for tier in ("hit", "opt")
            for name, uid in det2uid[tier].items()
        }
    )
    lh5.write(detector_uids, "detector_uids", evt_file, wo_mode="append")

    # forward the number of simulated primary events that remage stores at the
    # root of the stp file, so it can be summed across jobs at the cvt tier and
    # used to normalise the pdf histograms.
    lh5.write(
        lh5.read("number_of_simulated_events", stp_file),
        "number_of_simulated_events",
        evt_file,
        wo_mode="append",
    )

    with perf_block("move_to_cfs()"):
        move2cfs()

    print_stats()


if __name__ == "__main__":
    main()
