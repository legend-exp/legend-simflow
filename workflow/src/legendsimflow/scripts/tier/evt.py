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
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import numpy as np
from dbetto import AttrsDict
from dbetto.utils import load_dict
from lgdo import Array, Table, VectorOfVectors, lh5

from legendsimflow import nersc, patterns, spms_pars, utils
from legendsimflow import reboost as reboost_utils
from legendsimflow.awkward import ak_isin
from legendsimflow.metadata import encode_psd_usability, encode_usability
from legendsimflow.profile import make_profiler
from legendsimflow.tcm import merge_stp_n_opt_tcms_to_lh5

GEDS_ENERGY_THR_KEV = 25
SPMS_ENERGY_THR_PE = 0
BUFFER_LEN = "50*MB"
OFF = encode_usability("off")
ON = encode_usability("on")
VALID_PSD = encode_psd_usability("valid")

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

wildcards = args.wildcards
jobid = wildcards.jobid
stp_file = patterns.output_simjob_filename(
    args.config, tier="stp", simid=wildcards.simid, jobid=jobid
)
hit_file = {
    "opt": args.input.opt_file,
    "hit": args.input.hit_file,
}
evt_file = args.output[0]
log_file = args.log[0]
metadata = args.config.metadata
simstat_part_file = args.input.simstat_part_file
add_random_coincidences = args.params.add_random_coincidences
l200data = args.config.paths.l200data
usabilities = AttrsDict(load_dict(args.input.detector_usabilities[0]))

evt_file, move2cfs = nersc.make_on_scratch(args.config, evt_file)

# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
perf_block, print_stats, print_stats_since_last = make_profiler()

log.info("merging hit and opt TCMs")
with perf_block("merge_tcms()"):
    scintillator_uid = next(
        uid
        for uid, name in reboost_utils.get_remage_detector_uids(stp_file).items()
        if name == "liquid_argon"
    )

    merge_stp_n_opt_tcms_to_lh5(
        stp_file,
        hit_file["opt"],
        evt_file,
        scintillator_uid=scintillator_uid,
        buffer_len=BUFFER_LEN,
    )


# test that the evt tcm has the same amount of rows as the stp tcm

if lh5.read_n_rows("tcm", stp_file) != lh5.read_n_rows("tcm", evt_file):
    msg = (
        "stp and evt tcm should have same number of rows not "
        f"stp={lh5.read_n_rows('tcm', stp_file)}, "
        f"evt={lh5.read_n_rows('tcm', evt_file)}, "
        f"hit={lh5.read_n_rows('tcm', hit_file['hit'])}, "
        f"opt={lh5.read_n_rows('tcm', hit_file['opt'])}"
    )

    raise ValueError(msg)

# get the mapping of detector name to uid
# NOTE: we check on disk because we are not sure which tables were processed in
# the hit tiers
det2uid = {}
for tier in ("opt", "hit"):
    det2uid[tier] = {
        name: uid
        for uid, name in reboost_utils.get_remage_detector_uids(
            hit_file[tier], lh5_table="hit"
        ).items()
    }
    msg = f"found mapping name -> uid ({tier} tier): {det2uid[tier]}"
    log.debug(msg)


# little helper to simplify the code below
# TODO: move/fix in reboost
def _read_hits(tcm_ak, tier, field):
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
                data_ch = lh5.read(f"hit/{tab_name}/{field}", hit_file[tier], idx=rows)

        units = data_ch.attrs.get("units", None)
        data_ch = data_ch.view_as("ak")

        data_flat.append(data_ch)

    tcm_rows_concat = np.concatenate(tcm_rows)
    data_flat_concat = ak.concatenate(data_flat)[np.argsort(tcm_rows_concat)]

    data_unflat = ak.unflatten(data_flat_concat, counts)

    if units is not None:
        return ak.with_parameter(data_unflat, "units", units)
    return data_unflat


partitions = load_dict(simstat_part_file)[f"job_{jobid}"]

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
    on_spms_uids = sorted(
        uid
        for det_name, uid in det2uid["opt"].items()
        if (usabilities[runid].get(det_name) or {}).get("usability", "on") != "off"
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
        buffer_len=BUFFER_LEN,
        h5py_open_mode="a",
    )
    for chunk in it:
        unified_tcm = chunk.view_as("ak")
        out_table = Table(size=len(unified_tcm))

        # split the unified TCM in two, one for each tier. in this way we will be
        # able to read data from each tier
        tcm = {}
        for tier in ("opt", "hit"):
            mask = ak_isin(unified_tcm.table_key, det2uid[tier].values())
            tcm[tier] = unified_tcm[mask]

        # trigger table
        # -------------
        out_table.add_field("trigger", Table(size=len(unified_tcm)))

        # global fields that are constant over the full events
        # let's take them from the hit tier
        for constant_field in ["run", "period", "evtid"]:
            data = _read_hits(tcm, "hit", constant_field)

            # sanity check
            assert len(data) == len(tcm["hit"])

            # replace the awkward missing values with NaN for LH5 compatibility
            data = ak.fill_none(ak.firsts(data, axis=-1), np.nan)
            out_table.add_field(f"trigger/{constant_field}", Array(data))

        timestamp = _read_hits(tcm, "hit", "t0")
        timestamp = ak.fill_none(ak.firsts(timestamp, axis=-1), np.nan)
        out_table.add_field(
            "trigger/timestamp", Array(timestamp, attrs={"units": "ns"})
        )

        # HPGe table
        # ----------
        out_table.add_field("geds", Table(size=len(unified_tcm)))

        # first read usability and energy
        usability = _read_hits(tcm, "hit", "usability")
        psd_usability = _read_hits(tcm, "hit", "psd_usability")
        energy = _read_hits(tcm, "hit", "energy")

        # we want to only store hits from events in ON and AC detectors and above
        # our energy threshold
        hitsel = (usability != OFF) & (energy > GEDS_ENERGY_THR_KEV)

        # we want to still be able to know which detectors are ON (and not AC)
        out_table.add_field(
            "geds/is_good_channel", VectorOfVectors(usability[hitsel] == ON)
        )
        out_table.add_field(
            "geds/energy", VectorOfVectors(energy[hitsel], attrs={"units": "keV"})
        )
        # NOTE: the energy sum does not include AC detectors
        out_table.add_field(
            "geds/energy_sum",
            Array(
                ak.sum(energy[hitsel & (usability == ON)], axis=-1),
                attrs={"units": "keV"},
            ),
        )

        # fields to identify detectors and lookup stuff in the lower tiers
        out_table.add_field("geds/rawid", VectorOfVectors(tcm["hit"].table_key[hitsel]))
        out_table.add_field(
            "geds/hit_idx", VectorOfVectors(tcm["hit"].row_in_table[hitsel])
        )

        # PSD subtable
        out_table.add_field("geds/psd", Table(size=len(unified_tcm)))
        out_table.add_field(
            "geds/psd/is_good", VectorOfVectors(psd_usability[hitsel] == VALID_PSD)
        )

        aoe = _read_hits(tcm, "hit", "aoe")
        out_table.add_field("geds/psd/aoe", VectorOfVectors(aoe[hitsel]))
        out_table.add_field("geds/psd/has_aoe", VectorOfVectors(~np.isnan(aoe[hitsel])))

        is_ss = _read_hits(tcm, "hit", "is_single_site")
        out_table.add_field("geds/psd/is_single_site", VectorOfVectors(is_ss[hitsel]))

        # compute multiplicity
        geds_multiplicity = ak.sum(hitsel, axis=-1)
        out_table.add_field("geds/multiplicity", Array(geds_multiplicity))

        # SiPM table
        # ----------
        out_table.add_field("spms", Table(size=len(unified_tcm)))

        # also here, we exclude the non usable channels. this is in line with what
        # done in the evt tier in pygama
        usability = _read_hits(tcm, "opt", "usability")
        energy = _read_hits(tcm, "opt", "energy")
        chansel = usability != OFF
        # we also discard all pulses with amplitude below threshold
        pesel = energy > SPMS_ENERGY_THR_PE

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
        out_table.add_field("spms/energy", VectorOfVectors(energy_sel))

        is_saturated = _read_hits(tcm, "opt", "is_saturated")
        is_saturated_sel = is_saturated[chansel]
        # fill in Falses for events with no LAr edep
        empty_is_saturated = ak.Array([[False for _ in on_spms_uids]] * n_events)
        is_saturated_sel = ak.where(is_empty_opt, empty_is_saturated, is_saturated_sel)
        out_table.add_field("spms/is_saturated", VectorOfVectors(is_saturated_sel))

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
            "spms/time", VectorOfVectors(time_sel, attrs={"units": "ns"})
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
            out_table.add_field("spms/rc_energy", VectorOfVectors(rc_chunk.npe))
            out_table.add_field(
                "spms/rc_time", VectorOfVectors(rc_chunk.t0, attrs={"units": "ns"})
            )

        # total amount of light per event
        energy_sum = ak.sum(ak.sum(energy[pesel][chansel], axis=-1), axis=-1)
        out_table.add_field("spms/energy_sum", Array(energy_sum))

        # how many channels saw some light
        spms_multiplicity = ak.sum(ak.any(chansel & pesel, axis=-1), axis=-1)
        out_table.add_field("spms/multiplicity", Array(spms_multiplicity))

        # coincidences table
        # ------------------
        out_table.add_field("coincident", Table(size=len(unified_tcm)))

        # is there a signal in the HPGe array?
        out_table.add_field("coincident/geds", Array(geds_multiplicity > 0))

        # is there a signal in the LAr instrumentation?
        lar_veto = (spms_multiplicity >= 4) | (energy_sum >= 4)
        out_table.add_field("coincident/spms", Array(lar_veto))

        # now write down
        with perf_block("write_chunk()"):
            lh5.write(out_table, "evt", evt_file, wo_mode=evt_wo_mode)
            evt_wo_mode = "append"

    print_stats_since_last()

with perf_block("move_to_cfs()"):
    move2cfs()


print_stats()
