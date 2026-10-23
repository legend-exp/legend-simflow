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
import functools
import importlib
import logging
from collections.abc import Callable, Mapping, Sequence
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
from numpy.typing import ArrayLike
from reboost.optmap.convolve import NumdetStats, OptmapForConvolve
from reboost.units import units_conv_ak
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import metadata as mutils
from legendsimflow import nersc, utils
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


def load_light_model(model: Mapping) -> Callable:
    """Import the function of a light model and bind its keyword arguments.

    `model` is the ``model`` block of a light source, see
    :func:`legendsimflow.metadata.get_opt_scintillators`.
    """
    module, _, name = model["function"].rpartition(".")
    return functools.partial(
        getattr(importlib.import_module(module), name), **model["kwargs"]
    )


def model_photoelectrons(
    model: Callable,
    xloc: VectorOfVectors | ak.Array,
    yloc: VectorOfVectors | ak.Array,
    zloc: VectorOfVectors | ak.Array,
    num_scint_ph: ak.Array,
    sipm: str,
    map_scaling: float = 1,
    max_pes_per_hit: int = -1,
) -> tuple[ak.Array, np.ndarray]:
    """Sample the photoelectrons of each energy deposition from a light model.

    The model replaces the optical map of a light source. It is called once
    per chunk as ``model(xloc, yloc, zloc, sipm=sipm)``, with flat arrays of
    the global coordinates (mm) of all steps, and `sipm` the name of the SiPM, or ``"all"`` for all SiPMs together. It
    returns the probability that a scintillation photon emitted in each step
    gives a photoelectron, at unit channel efficiency.

    As with the optical map, the number of photoelectrons of each step is
    drawn from a Poisson distribution of mean `num_scint_ph` times the
    probability times `map_scaling`, and only the first `max_pes_per_hit`
    photoelectrons of each hit are kept, in step order, if `max_pes_per_hit`
    is positive.

    Returns
    -------
    the number of photoelectrons of each step, with the shape of
    `num_scint_ph`, and
    the mean number of photoelectrons of each hit at unit channel efficiency.
    """
    steps = ak.to_numpy(ak.num(num_scint_ph))
    prob = np.asarray(
        model(
            *(
                ak.to_numpy(ak.flatten(units_conv_ak(a, "mm")))
                for a in (xloc, yloc, zloc)
            ),
            sipm=sipm,
        ),
        dtype=np.float64,
    )
    if prob.shape != (ak.sum(steps),) or not np.all((prob >= 0) & (prob <= 1)):
        msg = (
            "the light model must return one detection probability in [0, 1] "
            "per energy deposition"
        )
        raise ValueError(msg)
    mean = ak.to_numpy(ak.flatten(num_scint_ph)) * prob

    n_pes = np.random.default_rng().poisson(map_scaling * mean)
    if max_pes_per_hit > 0:
        cum = np.cumsum(n_pes)
        # photoelectrons of the earlier steps of the same hit
        before = (
            cum
            - n_pes
            - np.repeat(np.concatenate([[0], cum])[np.cumsum(steps) - steps], steps)
        )
        n_pes = np.minimum(n_pes, np.maximum(max_pes_per_hit - before, 0))

    return ak.unflatten(n_pes, steps), ak.to_numpy(
        ak.sum(ak.unflatten(mean, steps), axis=1)
    )


def group_by_window(
    values: ak.Array | ArrayLike, window_idx: ArrayLike, n_windows: int
) -> ak.Array:
    """Group rows of scintillator hits by the event window they belong to.

    Parameters
    ----------
    values
        one entry per scintillator hit, possibly from several tables.
    window_idx
        index of the window of each hit, in ``[0, n_windows)``.
    n_windows
        number of windows. Windows without hits get an empty list.

    Returns
    -------
    array with one list per window, holding the entries of its hits in the
    order of `values`. Lists of `values` gain one level of nesting.
    """
    window_idx = np.asarray(window_idx, dtype=np.int64)
    order = np.argsort(window_idx, kind="stable")
    counts = np.bincount(window_idx, minlength=n_windows)
    return ak.unflatten(ak.Array(values)[order], counts)


@snakemake_compatible(
    mapping={
        "stp_file": "input.stp_file",
        "optmaps": lambda snakemake: [
            f"{key.removeprefix('optmap_')}={path}"
            for key, path in snakemake.input.items()
            if key.startswith("optmap_")
        ],
        "geom_file": "input.geom",
        "simstat_part_file": "input.simstat_part_file",
        "usability_file": "input.usability",
        "jobid": "wildcards.jobid",
        "opt_file": "output[0]",
        "log_file": "log[0]",
        "optmap_per_sipm": "params.optmap_per_sipm",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(description="Build the opt tier.")
    parser.add_argument("--stp-file", required=True, help="input stp tier file")
    parser.add_argument(
        "--optmap",
        dest="optmaps",
        action="append",
        default=[],
        metavar="SOURCE=FILE",
        help=(
            "optical map of a light source, needed with light_source: optmap "
            "for each source without a light model"
        ),
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
    simstat_part_file = nersc.dvs_ro(config, args.simstat_part_file)
    usability_map = AttrsDict(load_dict(nersc.dvs_ro(config, args.usability_file)))

    opt_file, move2cfs = nersc.make_on_scratch(config, opt_file)

    tier_opt_settings = get_tier_settings(config, "opt")
    photoelectron_resolution_sigma = tier_opt_settings.photoelectron_resolution_sigma
    time_resolution_in_ns = tier_opt_settings.time_resolution_in_ns
    buffer_len = tier_opt_settings.buffer_len
    store_expected_pes = tier_opt_settings.store_expected_pes
    light_source = tier_opt_settings.light_source
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

    # setup logging
    log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "tier-opt", parser, args)
    perf_block, print_perf, print_perf_last = reboost.make_profiler()

    # load the geometry and retrieve registered sensitive volume tables
    geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
    scintillators = pygeomtools.detectors.get_all_senstables(
        geom, type_filter="scintillator"
    )
    sipms = pygeomtools.detectors.get_all_senstables(geom, type_filter="optical")

    if light_source == "optmap":
        sources = mutils.get_opt_scintillators(config)
        optmap_files = {}
        for item in args.optmaps:
            name, _, path = item.partition("=")
            optmap_files[name] = nersc.dvs_ro(config, path)
        mapped = sorted(n for n, s in sources.items() if s.model is None)
        if sorted(optmap_files) != mapped:
            msg = (
                f"one optical map per light source without a light model is "
                f"needed: sources {mapped}, maps given for {sorted(optmap_files)}"
            )
            raise SimflowConfigError(msg, settings_block)
        light_models = {
            name: load_light_model(s.model)
            for name, s in sources.items()
            if s.model is not None
        }

        # fail early and loudly: without a matching table the output file
        # would never be created. tracked photons do not need these tables
        for name, source in sources.items():
            if source.table not in scintillators:
                msg = (
                    f"scintillator table {source.table} of light source {name} not "
                    f"found in {gdml_file}. scintillator tables in the geometry: "
                    f"{sorted(scintillators)}"
                )
                raise SimflowConfigError(msg, settings_block)
        uid2source = {scintillators[s.table].uid: name for name, s in sources.items()}

        stp_tables = lh5.ls(stp_file, "stp/*")
        # chunk size in rows, the smallest of the scintillator tables
        rows_per_chunk = min(
            (
                LH5Iterator(
                    stp_file, f"stp/{s.table}", buffer_len=buffer_len
                ).buffer_len
                for s in sources.values()
                if f"stp/{s.table}" in stp_tables
            ),
            default=1,
        )

    def write_sipm_hits(
        pe_times: ak.Array,
        is_saturated: ArrayLike,
        expected_pes: ArrayLike | None,
        hits: Table,
        sipm: str,
        sipm_uid: int,
        runid: str,
        usability: str,
    ) -> None:
        """Apply the SiPM response to the photoelectrons and write them.

        `pe_times` holds the sorted photoelectron times of each output row, in
        nanoseconds since the start of the event. `hits` holds the ``evtid``
        and ``t0`` fields of the output rows.
        """
        with perf_block("photoelectron_resolution()"):
            # smear_photoelectrons() expects a FWHM, convert from sigma
            pe_amps = reboost.spms.smear_photoelectrons(
                pe_times,
                photoelectron_resolution_sigma * 2 * np.sqrt(2 * np.log(2)),
            )

        if time_resolution_in_ns > 0:
            with perf_block("cluster_photoelectrons()"):
                pe_times, pe_amps = reboost.spms.cluster_photoelectrons(
                    pe_times, pe_amps, time_resolution_in_ns
                )

        with perf_block("write_hit_table_chunk()"):
            # relative to the hit t0, subtracted in float64
            pe_times = pe_times - hits.t0.view_as("ak")
            hits.add_field(
                "dt",
                VectorOfVectors(
                    ak.values_astype(pe_times, np.float32), attrs={"units": "ns"}
                ),
            )
            hits.add_field(
                "energy", VectorOfVectors(ak.values_astype(pe_amps, np.float32))
            )
            hits.add_field("is_saturated", Array(np.asarray(is_saturated)))
            if expected_pes is not None:
                hits.add_field(
                    "expected_pes",
                    Array(np.asarray(expected_pes, dtype=np.float32)),
                )

            _, period, run, _ = mutils.parse_runid(runid)
            field_vals = [period, run, mutils.encode_usability(usability)]
            for i, field in enumerate(["period", "run", "usability"]):
                hits.add_field(
                    field,
                    Array(np.full(shape=len(hits), fill_value=field_vals[i])),
                )

            reboost.write_hit_table_chunk(
                hits,
                "hit/" + ("spms" if sipm == "all" else sipm),
                opt_file,
                uid=sipm_uid,
            )

    def process_tracked_photons(
        iterator: Sequence[Table] | LH5Iterator,
        sipm: str,
        sipm_uid: int,
        runid: str,
        usability: str,
    ) -> None:
        """Write the photoelectrons of the photons tracked to one SiPM.

        Each row of the SiPM table gives one output row.
        """
        for lgdo_chunk in iterator:
            chunk = lgdo_chunk.view_as("ak")

            # the SiPM efficiency is applied while tracking: every recorded
            # photon is a photoelectron
            pe_times = ak.sort(chunk.time, axis=-1)
            if max_pes_per_hit > 0:
                # as with the map: saturated once the cap is reached
                is_saturated = ak.to_numpy(ak.num(pe_times) >= max_pes_per_hit)
                pe_times = pe_times[:, :max_pes_per_hit]
            else:
                is_saturated = np.full(len(chunk), fill_value=False, dtype=np.bool_)

            write_sipm_hits(
                pe_times,
                is_saturated,
                None,
                reboost.init_hit_table(lgdo_chunk),
                sipm,
                sipm_uid,
                runid,
                usability,
            )

    def process_optmap(
        windows: Mapping[str, tuple[int, np.ndarray]],
        n_windows: int,
        models: Mapping[str, OptmapForConvolve | Callable],
        sipm: str,
        sipm_uid: int,
        runid: str,
        usability: str,
    ) -> None:
        """Write the photoelectrons of the scintillation light seen by one SiPM.

        `windows` maps each light source to the first row of its table in
        this partition and the output row (window) of each following row. The
        light of all hits of a window is pooled before the SiPM response is
        applied. `models` holds the optical map or the light model of each
        source.
        """
        scalings = {}
        for name, source in sources.items():
            # constant for this SiPM, so resolve it once instead of per chunk
            scalings[name] = resolve_map_scaling(source.optmap_scaling_factor, sipm)
            msg = f"using optical map scaling factor {scalings[name]} for {name} and {sipm}"
            log.debug(msg)

        stats = {
            name: NumdetStats()
            for name, model in models.items()
            if isinstance(model, OptmapForConvolve)
        }

        # hits of all sources per window, as cumulative counts
        cum_rows = np.cumsum(
            np.bincount(
                np.concatenate([w for _, w in windows.values()]), minlength=n_windows
            )
        )

        w_start = 0
        while w_start < n_windows:
            # windows are never split, so a chunk can exceed rows_per_chunk
            base = cum_rows[w_start - 1] if w_start > 0 else 0
            w_end = max(
                w_start + 1,
                int(np.searchsorted(cum_rows, base + rows_per_chunk, side="right")),
            )
            w_end = min(w_end, n_windows)

            pe_times, expected_pes, t0, evtid, w_idx = ([] for _ in range(5))
            for source, (first_row, w) in windows.items():
                i, j = np.searchsorted(w, [w_start, w_end])
                if i == j:
                    continue

                material = sources[source].material
                lgdo_chunk = lh5.read(
                    f"stp/{sources[source].table}",
                    stp_file,
                    start_row=first_row + i,
                    n_rows=j - i,
                )
                chunk = lgdo_chunk.view_as("ak")

                with perf_block("emitted_scintillation_photons()"):
                    # TODO: we execute this multiple times...
                    scint_ph = reboost.spms.emitted_scintillation_photons(
                        chunk.edep, chunk.particle, material
                    )

                if source in stats:
                    with perf_block("number_of_detected_photoelectrons()"):
                        # return_stats=True also silences the per-chunk warnings
                        # of reboost about steps outside the map
                        *_output, _stats = (
                            reboost.spms.number_of_detected_photoelectrons(
                                chunk.xloc,
                                chunk.yloc,
                                chunk.zloc,
                                scint_ph,
                                models[source],
                                sipm,
                                map_scaling=scalings[source],
                                max_pes_per_hit=max_pes_per_hit,
                                return_pes_expectation_value=store_expected_pes,
                                return_stats=True,
                            )
                        )
                    stats[source] += _stats

                    # reboost appends the expectation after the other outputs
                    if store_expected_pes:
                        *_output, _expected = _output
                        expected_pes.append(np.asarray(_expected, dtype=np.float64))
                    nr_pe = _output[0]
                else:
                    with perf_block("model_photoelectrons()"):
                        # the LGDO fields carry the units
                        nr_pe, _expected = model_photoelectrons(
                            models[source],
                            lgdo_chunk.xloc,
                            lgdo_chunk.yloc,
                            lgdo_chunk.zloc,
                            scint_ph,
                            sipm,
                            map_scaling=scalings[source],
                            max_pes_per_hit=max_pes_per_hit,
                        )
                    if store_expected_pes:
                        expected_pes.append(_expected)
                # a saturated hit has exactly max_pes_per_hit photoelectrons, so
                # its window is saturated by the check on the pooled light below

                with perf_block("photoelectron_times()"):
                    pe_times.append(
                        ak.without_parameters(
                            reboost.spms.photoelectron_times(
                                nr_pe, chunk.particle, chunk.time, material
                            )
                        )
                    )

                hits = reboost.init_hit_table(lgdo_chunk)
                t0.append(hits.t0.nda)
                evtid.append(hits.evtid.nda)
                w_idx.append(w[i:j] - w_start)

            n = w_end - w_start
            w_idx = np.concatenate(w_idx)

            with perf_block("pool_photoelectrons()"):
                # photoelectron_times() does not guarantee time ordering
                pooled = ak.sort(
                    ak.flatten(
                        group_by_window(ak.concatenate(pe_times), w_idx, n), axis=2
                    ),
                    axis=-1,
                )
                if max_pes_per_hit > 0:
                    saturated = ak.to_numpy(ak.num(pooled) >= max_pes_per_hit)
                    pooled = pooled[:, :max_pes_per_hit]
                else:
                    saturated = np.full(n, fill_value=False, dtype=np.bool_)

                expected = None
                if store_expected_pes:
                    expected = ak.to_numpy(
                        ak.sum(
                            group_by_window(np.concatenate(expected_pes), w_idx, n),
                            axis=1,
                        )
                    )

                window_hits = Table(size=n)
                window_hits.add_field(
                    "t0",
                    Array(
                        ak.to_numpy(
                            ak.min(
                                group_by_window(np.concatenate(t0), w_idx, n),
                                axis=1,
                                mask_identity=False,
                            )
                        ),
                        attrs={"units": "ns"},
                    ),
                )
                window_hits.add_field(
                    "evtid",
                    Array(
                        # every window has at least one hit
                        ak.to_numpy(
                            group_by_window(np.concatenate(evtid), w_idx, n)[:, 0]
                        )
                    ),
                )

            write_sipm_hits(
                pooled,
                saturated,
                expected,
                window_hits,
                sipm,
                sipm_uid,
                runid,
                usability,
            )
            w_start = w_end

        for name, s in stats.items():
            tot = s.energy_looped
            for counts, where in (
                (s.energy_oob + s.energy_no_stats, "outside the map"),
                (s.energy_zero, "in bins with zero detection probability"),
            ):
                pct = 100 * counts / tot if tot > 0 else 0.0
                msg = (
                    f"optical map {name}/{sipm}: {pct:.2f}% of the energy deposited "
                    f"in {name} is {where}"
                )
                log.log(logging.WARNING if pct > 1 else logging.DEBUG, msg)

    def get_usability(runid: str, sipm: str) -> str:
        usability = usability_map[runid].get(sipm)
        if usability is None:
            msg = f"usability not found for {sipm} in {runid}, defaulting to on"
            log.warning(msg)
            usability = "on"
        return usability

    partitions = load_dict(simstat_part_file)[f"job_{jobid}"]

    # load TCM, to be used to chunk the event statistics according to the run partitioning
    msg = "loading TCM"
    log.debug(msg)
    tcm = lh5.read_as("tcm", stp_file, library="ak")

    # in combined mode there is a single map per source, so pre-load them once
    # for a little speed up. in per-SiPM mode the maps are loaded on demand
    if light_source == "optmap" and not optmap_per_sipm:
        optmaps_all = light_models | {
            name: reboost.spms.load_optmap(f, "all") for name, f in optmap_files.items()
        }

    # loop over the partitions for this file
    for runid_idx, (runid, evt_idx_range) in enumerate(partitions.items()):
        msg = (
            f"processing partition corresponding to {runid} "
            f"[{runid_idx + 1}/{len(partitions)}], event range {evt_idx_range}"
        )
        log.info(msg)

        if light_source == "tracked_photons":
            # as in the hit tier, each row of a SiPM table gives one output row
            for sipm in sorted(sipms):
                sipm_uid = sipms[sipm].uid
                i_start, n_entries = reboost.get_rows_in_event_range(
                    tcm, sipm_uid, *evt_idx_range
                )
                if n_entries == 0:
                    continue

                msg = f"processing the tracked photons of SiPM {sipm}"
                log.debug(msg)

                process_tracked_photons(
                    LH5Iterator(
                        stp_file,
                        f"stp/{sipm}",
                        i_start=i_start,
                        n_entries=n_entries,
                        buffer_len=buffer_len,
                        field_mask=["evtid", "t0", "time"],
                    ),
                    sipm,
                    sipm_uid,
                    runid,
                    get_usability(runid, sipm),
                )
                print_perf_last()
            continue

        # each event window (TCM row) with a scintillator hit gives one output
        # row, where the light of all its scintillator hits is pooled
        first, last = evt_idx_range
        tcm_part = tcm[first : last + 1]
        keys = ak.to_numpy(ak.flatten(tcm_part.table_key))
        rows = ak.to_numpy(ak.flatten(tcm_part.row_in_table))
        in_window = np.repeat(
            np.arange(len(tcm_part)), ak.to_numpy(ak.num(tcm_part.table_key))
        )

        is_scint = np.isin(keys, list(uid2source))
        has_scint = np.zeros(len(tcm_part), dtype=np.bool_)
        has_scint[in_window[is_scint]] = True
        out_row = np.cumsum(has_scint) - 1
        n_windows = int(has_scint.sum())

        windows = {}
        for uid, source in uid2source.items():
            sel = keys == uid
            if not np.any(sel):
                continue
            # the TCM orders the hits of a window by time: sort them by row,
            # so that each chunk is read as one block of rows
            order = np.argsort(rows[sel], kind="stable")
            src_rows = rows[sel][order]
            src_w = out_row[in_window[sel]][order]
            if np.any(np.diff(src_rows) != 1) or np.any(np.diff(src_w) < 0):
                msg = (
                    f"the rows of stp/{sources[source].table} in event range "
                    f"{evt_idx_range} are not one block in event order"
                )
                raise RuntimeError(msg)
            windows[source] = (int(src_rows[0]), src_w)

        msg = (
            f"found {n_windows} event windows with scintillator hits, "
            f"{ {n: len(w) for n, (_, w) in windows.items()} } hits per light source"
        )
        log.info(msg)
        if n_windows == 0:
            continue

        if optmap_per_sipm:
            for sipm in sorted(sipms):
                msg = f"applying optical maps for SiPM {sipm}"
                log.debug(msg)

                with perf_block("load_optmap()"):
                    # the per-channel maps are too large to keep all of them
                    # resident for the whole job: load and release them here
                    optmaps = light_models | {
                        name: reboost.spms.load_optmap(f, sipm)
                        for name, f in optmap_files.items()
                    }

                process_optmap(
                    windows,
                    n_windows,
                    optmaps,
                    sipm,
                    sipms[sipm].uid,
                    runid,
                    get_usability(runid, sipm),
                )
                del optmaps
                print_perf_last()
        else:
            log.debug("applying sum optical maps and light models")

            process_optmap(
                windows,
                n_windows,
                optmaps_all,
                "all",
                # the smallest UID of the light sources, independent of their
                # order in the configuration
                min(uid2source),
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
        for sipm in sorted(sipms):
            if f"hit/{sipm}" in tables:
                continue
            process_tracked_photons([no_photons], sipm, sipms[sipm].uid, runid, "on")

    log.debug("building the TCM")
    reboost.build_remage_tcm(opt_file, opt_file)

    with perf_block("move_to_cfs()"):
        move2cfs()

    print_perf()


if __name__ == "__main__":
    main()
