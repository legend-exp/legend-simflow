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

"""Build the photoelectron tables of the opt tier.

The opt tier holds, for each SiPM channel, the detected photoelectrons.
Photoelectrons come either from the convolution of the energy depositions in
the scintillators (the light sources) with their optical maps or light models
(``light_source: optmap``), or from the optical photons tracked by Geant4 up to the SiPMs
(``light_source: tracked_photons``). With optical maps, the light of all
scintillator hits in the same event window (row of the stp TCM) is pooled and
gives one output row. With tracked photons, each hit of a SiPM gives one output
row. Times are in nanoseconds, amplitudes in photoelectrons.
"""

from __future__ import annotations

import functools
import importlib
import inspect
import logging
from collections.abc import Callable, Iterable, Mapping
from contextlib import AbstractContextManager, nullcontext
from pathlib import Path
from typing import TYPE_CHECKING

import awkward as ak
import lh5
import numpy as np
import reboost
import reboost.spms
from lgdo import Array, Table, VectorOfVectors
from numpy.typing import ArrayLike, NDArray
from reboost.optmap.convolve import NumdetStats, OptmapForConvolve
from reboost.units import units_conv_ak

from legendsimflow import metadata as mutils

if TYPE_CHECKING:
    from pyg4ometry.geant4 import Registry

log = logging.getLogger(__name__)


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


def get_usability(
    usability_map: Mapping[str, Mapping[str, str]], runid: str, channel: str
) -> str:
    """Usability of `channel` in `runid`, ``"on"`` with a warning if missing."""
    usability = usability_map[runid].get(channel)
    if usability is None:
        msg = f"usability not found for {channel} in {runid}, defaulting to on"
        log.warning(msg)
        usability = "on"
    return usability


def tracked_photons_to_pes(
    chunk: ak.Array, max_pes_per_hit: int
) -> tuple[ak.Array, NDArray[np.bool_]]:
    """Photoelectron times of the optical photons tracked up to a SiPM.

    The SiPM efficiency is applied while tracking, so every recorded photon is
    a photoelectron.

    Parameters
    ----------
    chunk
        Rows of a SiPM step table, with the arrival times of the photons in
        the ``time`` field (one list per hit, in nanoseconds).
    max_pes_per_hit
        Keep at most this many photoelectrons per hit, the earliest ones. A
        hit that reaches the cap is flagged as saturated. Zero or negative
        keeps all photoelectrons.

    Returns
    -------
    pe_times
        The photoelectron times of each hit, in increasing order.
    is_saturated
        One flag per hit.
    """
    pe_times = ak.sort(chunk.time, axis=-1)
    if max_pes_per_hit > 0:
        is_saturated = ak.to_numpy(ak.num(pe_times) >= max_pes_per_hit)
        pe_times = pe_times[:, :max_pes_per_hit]
    else:
        is_saturated = np.full(len(chunk), fill_value=False, dtype=np.bool_)
    return pe_times, is_saturated


def optmap_to_pes(
    chunk: ak.Array,
    optmap: OptmapForConvolve,
    sipm: str,
    *,
    material: str,
    map_scaling: float,
    max_pes_per_hit: int,
    store_expected_pes: bool,
    perf_block: Callable[[str], AbstractContextManager] = nullcontext,
) -> tuple[ak.Array, NDArray[np.bool_], NDArray[np.floating] | None, NumdetStats]:
    """Sample the photoelectrons of the energy depositions in a scintillator.

    The number of scintillation photons emitted by each step follows from its
    energy deposition, particle type and the light yield of `material`. The
    optical map gives the probability that a photon emitted at the step
    position reaches `sipm`. The emission times follow the scintillation time
    profile of `material`.

    Parameters
    ----------
    chunk
        Rows of a scintillator step table (fields ``edep``, ``particle``,
        ``xloc``, ``yloc``, ``zloc`` and ``time``, one list per hit).
    optmap
        The optical map, already loaded with :func:`reboost.spms.load_optmap`.
    sipm
        Name of the SiPM, or ``"all"`` for the combined map.
    material
        Scintillation material, as known to reboost: ``lar``, ``pen`` or
        ``fiber``.
    map_scaling
        Factor applied to the detection probabilities of the map.
    max_pes_per_hit
        Cap on the photoelectrons sampled per hit. A hit that reaches the cap
        is flagged as saturated. Zero or negative disables the cap.
    store_expected_pes
        Also return the expected number of photoelectrons of each hit.
    perf_block
        Context manager of a :func:`reboost.profile.make_profiler`, used to time the
        reboost calls.

    Returns
    -------
    pe_times
        The photoelectron times of each hit, in increasing order.
    is_saturated
        One flag per hit.
    expected_pes
        The expected number of photoelectrons of each hit, or ``None``.
    stats
        Energy depositions outside the map or in bins without detection
        probability, summed over the chunk.
    """
    with perf_block("emitted_scintillation_photons()"):
        scint_ph = reboost.spms.emitted_scintillation_photons(
            chunk.edep, chunk.particle, material
        )

    with perf_block("number_of_detected_photoelectrons()"):
        # return_stats=True also silences the per-chunk warnings of reboost
        # about steps outside the map
        *output, stats = reboost.spms.number_of_detected_photoelectrons(
            chunk.xloc,
            chunk.yloc,
            chunk.zloc,
            scint_ph,
            optmap,
            sipm,
            map_scaling=map_scaling,
            max_pes_per_hit=max_pes_per_hit,
            return_pes_expectation_value=store_expected_pes,
            return_stats=True,
        )

    # reboost appends the expectation after the other outputs
    if store_expected_pes:
        *output, expected_pes = output
    else:
        expected_pes = None
    if max_pes_per_hit > 0:
        nr_pe, is_saturated = output
    else:
        (nr_pe,) = output
        is_saturated = np.full(len(chunk), fill_value=False, dtype=np.bool_)

    with perf_block("photoelectron_times()"):
        pe_times = reboost.spms.photoelectron_times(
            nr_pe, chunk.particle, chunk.time, material
        )
        # the photoelectron_times() processor does not guarantee time ordering
        pe_times = ak.sort(pe_times, axis=-1)

    return pe_times, is_saturated, expected_pes, stats


def load_light_model(model: Mapping, geometry: Registry | None = None) -> Callable:
    """Import the function of a light model and bind its keyword arguments.

    `model` is the ``model`` block of a light source, see
    :func:`.metadata.get_opt_scintillators`. A function with a ``geometry``
    keyword argument also gets `geometry`, the simulated GDML geometry.
    """
    module, _, name = model["function"].rpartition(".")
    func = getattr(importlib.import_module(module), name)
    kwargs = dict(model["kwargs"])
    if "geometry" in inspect.signature(func).parameters:
        if geometry is None:
            msg = f"the light model {model['function']} needs the GDML geometry"
            raise ValueError(msg)
        kwargs["geometry"] = geometry
    return functools.partial(func, **kwargs)


def load_light_sources(
    stp_file: str | Path,
    sources: Mapping[str, Mapping],
    optmap_files: Mapping[str, str | Path],
    geometry: Registry,
    *,
    buffer_len: str | int,
    preload_combined: bool,
) -> tuple[dict[str, str | Path | OptmapForConvolve | Callable], int]:
    """Prepare the optical maps and light models of the light sources.

    A source without a table in `stp_file` has no hits: it gets a warning, and
    its map is not loaded.

    Parameters
    ----------
    stp_file
        The stp tier file.
    sources
        The light sources, see :func:`.metadata.get_opt_scintillators`.
    optmap_files
        The optical map file of each source without a light model.
    geometry
        The simulated GDML geometry, for the light models that need it.
    buffer_len
        The LH5 read chunk size of the ``opt`` tier settings, e.g.
        ``"500*MB"``.
    preload_combined
        Load the map of all SiPMs together (``all``) of each source with hits
        here. Otherwise the maps are loaded by
        :func:`write_scintillation_photoelectrons`, one SiPM at a time.

    Returns
    -------
    optmaps
        For each source, its optical map or map file, or its light model.
    rows_per_chunk
        Number of scintillator hits to read at once: the smallest that
        `buffer_len` gives for the tables of the sources with hits.
    """
    stp_tables = lh5.ls(stp_file, "stp/*")
    with_hits = []
    for name, source in sources.items():
        if f"stp/{source['table']}" in stp_tables:
            with_hits.append(name)
        else:
            msg = (
                f"table {source['table']} of light source {name} not found in "
                f"{stp_file}. possibly because it was not read-out or there "
                "were no hits recorded"
            )
            log.warning(msg)

    optmaps = {}
    for name, source in sources.items():
        if source["model"] is not None:
            optmaps[name] = load_light_model(source["model"], geometry)
        elif preload_combined and name in with_hits:
            optmaps[name] = reboost.spms.load_optmap(optmap_files[name], "all")
        else:
            optmaps[name] = optmap_files[name]

    rows_per_chunk = min(
        (
            lh5.LH5Iterator(
                str(stp_file), f"stp/{sources[name]['table']}", buffer_len=buffer_len
            ).buffer_len
            for name in with_hits
        ),
        default=1,
    )
    return optmaps, rows_per_chunk


def model_photoelectrons(
    model: Callable,
    xloc: VectorOfVectors | ak.Array,
    yloc: VectorOfVectors | ak.Array,
    zloc: VectorOfVectors | ak.Array,
    num_scint_ph: ak.Array,
    map_scaling: float = 1,
    max_pes_per_hit: int = -1,
) -> tuple[ak.Array, NDArray[np.floating]]:
    """Sample the photoelectrons of each energy deposition from a light model.

    The model replaces the optical map of a light source. It is called once
    per chunk as ``model(xloc, yloc, zloc)``, with flat arrays of the global
    coordinates (mm) of all steps. It returns the probability that a
    scintillation photon emitted in each step gives a photoelectron in any
    SiPM, at unit channel efficiency.

    `xloc`, `yloc` and `zloc` are converted to mm from the units in their
    attributes (LGDO) or parameters (awkward array). Positions without units
    are taken in mm: pass the step table columns, not their default awkward
    view, which drops the units (metres in remage).

    As with the optical map, the number of photoelectrons of each step is
    drawn from a Poisson distribution of mean `num_scint_ph` times the
    probability times `map_scaling`, and only the first `max_pes_per_hit`
    photoelectrons of each hit are kept, in step order, if `max_pes_per_hit`
    is positive.

    Returns
    -------
    n_pes
        The number of photoelectrons of each step, with the shape of
        `num_scint_ph`.
    expected_pes
        The mean number of photoelectrons of each hit at unit channel
        efficiency.
    """
    steps = ak.to_numpy(ak.num(num_scint_ph))
    prob = np.asarray(
        model(
            *(
                ak.to_numpy(ak.flatten(units_conv_ak(a, "mm")))
                for a in (xloc, yloc, zloc)
            ),
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


def light_model_to_pes(
    stp_chunk: Table,
    model: Callable,
    *,
    material: str,
    map_scaling: float,
    max_pes_per_hit: int,
    perf_block: Callable[[str], AbstractContextManager] = nullcontext,
) -> tuple[ak.Array, NDArray[np.floating]]:
    """Sample the photoelectrons of the energy depositions with a light model.

    As :func:`optmap_to_pes`, with the detection probabilities given by
    `model` (see :func:`model_photoelectrons`) instead of an optical map.
    `stp_chunk` is a scintillator step table, whose position fields carry
    their units.

    Returns the photoelectron times of each hit, in increasing order, and the
    expected number of photoelectrons of each hit.
    """
    chunk = stp_chunk.view_as("ak")
    with perf_block("emitted_scintillation_photons()"):
        scint_ph = reboost.spms.emitted_scintillation_photons(
            chunk.edep, chunk.particle, material
        )

    with perf_block("model_photoelectrons()"):
        nr_pe, expected_pes = model_photoelectrons(
            model,
            stp_chunk.xloc,
            stp_chunk.yloc,
            stp_chunk.zloc,
            scint_ph,
            map_scaling=map_scaling,
            max_pes_per_hit=max_pes_per_hit,
        )

    with perf_block("photoelectron_times()"):
        pe_times = reboost.spms.photoelectron_times(
            nr_pe, chunk.particle, chunk.time, material
        )
        pe_times = ak.sort(pe_times, axis=-1)

    return pe_times, expected_pes


def apply_photoelectron_resolution(
    pe_times: ak.Array,
    sigma_in_pe: float,
    time_resolution_in_ns: float,
    perf_block: Callable[[str], AbstractContextManager] = nullcontext,
) -> tuple[ak.Array, ak.Array]:
    """Smear the photoelectron amplitudes and merge the unresolved ones.

    Each photoelectron gets an amplitude drawn from a Gaussian of mean one
    and standard deviation `sigma_in_pe` (in photoelectrons). Photoelectrons
    closer in time than `time_resolution_in_ns` are then merged into one, at
    the time of the first and with the summed amplitude. A time resolution of
    zero or less skips the merging.

    Returns the times and the amplitudes, with one list per hit.
    """
    with perf_block("photoelectron_resolution()"):
        # smear_photoelectrons() expects a FWHM, convert from sigma
        pe_amps = reboost.spms.smear_photoelectrons(
            pe_times, sigma_in_pe * 2 * np.sqrt(2 * np.log(2))
        )

    if time_resolution_in_ns > 0:
        with perf_block("cluster_photoelectrons()"):
            pe_times, pe_amps = reboost.spms.cluster_photoelectrons(
                pe_times, pe_amps, time_resolution_in_ns
            )

    return pe_times, pe_amps


def build_opt_table(
    stp_chunk: Table,
    pe_times: ak.Array,
    pe_amps: ak.Array,
    is_saturated: NDArray[np.bool_],
    expected_pes: NDArray[np.floating] | None,
    runid: str,
    usability: str,
) -> Table:
    """Assemble the opt tier table of one SiPM for a chunk of output rows.

    The table carries the event fields of `stp_chunk`, a step table or a table
    of event windows with ``evtid`` and ``t0`` (see
    :func:`reboost.init_hit_table`), and, per row: ``dt`` (photoelectron times
    relative to the hit time ``t0``, in nanoseconds), ``energy`` (the
    amplitudes, in photoelectrons), ``is_saturated``, ``expected_pes`` (only
    if given), and the ``period``, ``run`` and encoded ``usability`` of the
    data run the row belongs to.
    """
    out_table = reboost.init_hit_table(stp_chunk)

    # relative to the hit t0, subtracted in float64
    pe_times = pe_times - out_table.t0.view_as("ak")
    out_table.add_field(
        "dt",
        VectorOfVectors(ak.values_astype(pe_times, np.float32), attrs={"units": "ns"}),
    )
    out_table.add_field(
        "energy", VectorOfVectors(ak.values_astype(pe_amps, np.float32))
    )
    out_table.add_field("is_saturated", Array(is_saturated))
    if expected_pes is not None:
        out_table.add_field(
            "expected_pes", Array(np.asarray(expected_pes, dtype=np.float32))
        )

    _, period, run, _ = mutils.parse_runid(runid)
    field_vals = [period, run, mutils.encode_usability(usability)]
    for i, field in enumerate(["period", "run", "usability"]):
        out_table.add_field(
            field, Array(np.full(shape=len(out_table), fill_value=field_vals[i]))
        )

    return out_table


def log_optmap_coverage(stats: NumdetStats, sipm: str, volume_name: str) -> None:
    """Log the fraction of the energy deposited where the map gives no light.

    Energy outside the map or in bins with zero detection probability is
    reported as a warning above one percent, as a debug message otherwise.
    """
    tot = stats.energy_looped
    for counts, where in (
        (stats.energy_oob + stats.energy_no_stats, "outside the map"),
        (stats.energy_zero, "in bins with zero detection probability"),
    ):
        pct = 100 * counts / tot if tot > 0 else 0.0
        msg = (
            f"optical map {sipm}: {pct:.2f}% of the energy deposited "
            f"in {volume_name} is {where}"
        )
        log.log(logging.WARNING if pct > 1 else logging.DEBUG, msg)


def empty_photon_table() -> Table:
    """A SiPM step table with no rows, to give a table to a SiPM without photons."""
    return Table(
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


def group_by_window(
    values: ak.Array | ArrayLike, window_idx: ArrayLike, n_windows: int
) -> ak.Array:
    """Group rows of scintillator hits by the event window they belong to.

    Parameters
    ----------
    values
        One entry per scintillator hit, possibly from several tables.
    window_idx
        Index of the window of each hit, in ``[0, n_windows)``.
    n_windows
        Number of windows. Windows without hits get an empty list.

    Returns
    -------
    ak.Array
        One list per window, holding the entries of its hits in the order of
        `values`. Lists of `values` gain one level of nesting.
    """
    window_idx = np.asarray(window_idx, dtype=np.int64)
    order = np.argsort(window_idx, kind="stable")
    counts = np.bincount(window_idx, minlength=n_windows)
    return ak.unflatten(ak.Array(values)[order], counts)


def find_event_windows(
    tcm: ak.Array, uid2source: Mapping[int, str]
) -> tuple[dict[str, tuple[int, NDArray[np.int64]]], int]:
    """Find the event windows with light from the scintillators.

    Each event window (TCM row) with at least one hit of a light source gives
    one output row of the opt tier.

    Parameters
    ----------
    tcm
        The rows of the stp TCM to process (fields ``table_key`` and
        ``row_in_table``).
    uid2source
        The light sources, by the UID of their table.

    Returns
    -------
    windows
        For each light source with hits, the first row of its table in `tcm`
        and the output row (window) of each hit from there on.
    n_windows
        The number of event windows with scintillator hits.
    """
    keys = ak.to_numpy(ak.flatten(tcm.table_key))
    rows = ak.to_numpy(ak.flatten(tcm.row_in_table))
    in_window = np.repeat(np.arange(len(tcm)), ak.to_numpy(ak.num(tcm.table_key)))

    is_scint = np.isin(keys, list(uid2source))
    has_scint = np.zeros(len(tcm), dtype=np.bool_)
    has_scint[in_window[is_scint]] = True
    out_row = np.cumsum(has_scint) - 1

    windows = {}
    for uid, source in uid2source.items():
        sel = keys == uid
        if not np.any(sel):
            continue
        # the TCM orders the hits of a window by time: sort them by row, so
        # that each chunk is read as one block of rows
        order = np.argsort(rows[sel], kind="stable")
        src_rows = rows[sel][order]
        src_w = out_row[in_window[sel]][order]
        if np.any(np.diff(src_rows) != 1) or np.any(np.diff(src_w) < 0):
            msg = f"the rows of light source {source} are not one block in event order"
            raise RuntimeError(msg)
        windows[source] = (int(src_rows[0]), src_w)

    return windows, int(has_scint.sum())


def pool_photoelectrons(
    pe_times: ak.Array,
    expected_pes: NDArray[np.floating] | None,
    t0: ArrayLike,
    evtid: ArrayLike,
    window_idx: ArrayLike,
    n_windows: int,
    max_pes_per_hit: int,
) -> tuple[Table, ak.Array, NDArray[np.bool_], NDArray[np.floating] | None]:
    """Pool the photoelectrons of the scintillator hits of each event window.

    Parameters
    ----------
    pe_times
        The photoelectron times of each hit, in nanoseconds.
    expected_pes
        The expected number of photoelectrons of each hit, or ``None``.
    t0
        The time of each hit, in nanoseconds.
    evtid
        The Geant4 event of each hit.
    window_idx
        Index of the window of each hit, in ``[0, n_windows)``. Every window
        has at least one hit.
    n_windows
        Number of windows.
    max_pes_per_hit
        Keep at most this many photoelectrons per window, the earliest ones.
        A window that reaches the cap is flagged as saturated. Zero or
        negative keeps all photoelectrons.

    Returns
    -------
    windows
        The ``evtid`` and ``t0`` (the earliest hit) of each window.
    pe_times
        The photoelectron times of each window, in increasing order.
    is_saturated
        One flag per window.
    expected_pes
        The expected number of photoelectrons of each window, or ``None``.
    """
    pooled = ak.sort(
        ak.flatten(group_by_window(pe_times, window_idx, n_windows), axis=2), axis=-1
    )
    if max_pes_per_hit > 0:
        is_saturated = ak.to_numpy(ak.num(pooled) >= max_pes_per_hit)
        pooled = pooled[:, :max_pes_per_hit]
    else:
        is_saturated = np.full(n_windows, fill_value=False, dtype=np.bool_)

    if expected_pes is not None:
        expected_pes = ak.to_numpy(
            ak.sum(group_by_window(expected_pes, window_idx, n_windows), axis=1)
        )

    windows = Table(size=n_windows)
    t0 = ak.min(group_by_window(t0, window_idx, n_windows), axis=1, mask_identity=False)
    windows.add_field("t0", Array(ak.to_numpy(t0), attrs={"units": "ns"}))
    evtid = group_by_window(evtid, window_idx, n_windows)[:, 0]
    windows.add_field("evtid", Array(ak.to_numpy(evtid)))

    return windows, pooled, is_saturated, expected_pes


def _write_chunk(
    stp_chunk: Table,
    pe_times: ak.Array,
    is_saturated: NDArray[np.bool_],
    expected_pes: NDArray[np.floating] | None,
    out_file: str | Path,
    sipm: str,
    sipm_uid: int,
    runid: str,
    usability: str,
    settings: Mapping,
    perf_block: Callable[[str], AbstractContextManager],
) -> None:
    """Apply the SiPM response to the photoelectrons and append them to `out_file`."""
    pe_times, pe_amps = apply_photoelectron_resolution(
        pe_times,
        settings["photoelectron_resolution_sigma"],
        settings["time_resolution_in_ns"],
        perf_block=perf_block,
    )

    with perf_block("write_hit_table_chunk()"):
        out_table = build_opt_table(
            stp_chunk, pe_times, pe_amps, is_saturated, expected_pes, runid, usability
        )
        reboost.write_hit_table_chunk(
            out_table,
            "hit/" + ("spms" if sipm == "all" else sipm),
            out_file,
            uid=sipm_uid,
        )


def write_tracked_photoelectrons(
    chunks: Iterable[Table],
    out_file: str | Path,
    sipm: str,
    sipm_uid: int,
    runid: str,
    usability: str,
    *,
    settings: Mapping,
    perf_block: Callable[[str], AbstractContextManager] = nullcontext,
) -> None:
    """Write the photoelectrons of the photons tracked up to one SiPM.

    Each row of `chunks`, the step table of the SiPM, gives one output row,
    appended to ``hit/<sipm>`` in `out_file`.

    Parameters
    ----------
    chunks
        The step table rows to process, as returned by
        :class:`~lh5.io.iterator.LH5Iterator`.
    out_file
        The opt tier file.
    sipm
        Name of the SiPM.
    sipm_uid
        UID of the SiPM.
    runid
        The data run the rows belong to, as in the simulation statistics
        partition file.
    usability
        Usability of the SiPM in `runid`.
    settings
        The ``opt`` tier settings (see :func:`.metadata.get_tier_settings`):
        ``max_pes_per_hit_per_sipm``, ``photoelectron_resolution_sigma`` and
        ``time_resolution_in_ns``.
    perf_block
        Context manager of a :func:`reboost.profile.make_profiler`, used to time the
        reboost calls.
    """
    for stp_chunk in chunks:
        pe_times, is_saturated = tracked_photons_to_pes(
            stp_chunk.view_as("ak"), settings["max_pes_per_hit_per_sipm"]
        )
        _write_chunk(
            stp_chunk,
            pe_times,
            is_saturated,
            None,
            out_file,
            sipm,
            sipm_uid,
            runid,
            usability,
            settings,
            perf_block,
        )


def write_scintillation_photoelectrons(
    stp_file: str | Path,
    windows: Mapping[str, tuple[int, NDArray[np.int64]]],
    n_windows: int,
    sources: Mapping[str, Mapping],
    optmaps: Mapping[str, str | Path | OptmapForConvolve | Callable],
    out_file: str | Path,
    sipm: str,
    sipm_uid: int,
    runid: str,
    usability: str,
    *,
    settings: Mapping,
    rows_per_chunk: int,
    perf_block: Callable[[str], AbstractContextManager] = nullcontext,
) -> dict[str, NumdetStats]:
    """Write the photoelectrons of the scintillation light seen by one SiPM.

    The light of all scintillator hits in an event window is pooled before
    the SiPM response is applied: the photoelectron cap, the amplitude
    smearing and the time clustering act on the combined light. Each window
    gives one output row, appended to ``hit/<sipm>`` in `out_file`
    (``hit/spms`` for the combined map).

    Parameters
    ----------
    stp_file
        The stp tier file.
    windows
        For each light source with hits, the first row of its table and the
        output row of each hit from there on, see :func:`find_event_windows`.
    n_windows
        Number of output rows.
    sources
        The light sources (``table``, ``material`` and
        ``optmap_scaling_factor``), see
        :func:`.metadata.get_opt_scintillators`.
    optmaps
        The optical map of each light source, its file name, or the light
        model that replaces it (see :func:`load_light_model`). Maps given by
        file name are loaded here, only for the sources with hits in
        `windows`, and released on return (the per-channel maps are too large
        to keep all of them in memory).
    out_file
        The opt tier file.
    sipm
        Name of the SiPM, or ``"all"`` for the combined map.
    sipm_uid
        UID of the output table, the one of the SiPM or, for the combined map,
        of a light source.
    runid
        The data run the windows belong to, as in the simulation statistics
        partition file.
    usability
        Usability of the SiPM in `runid`.
    settings
        The ``opt`` tier settings (see :func:`.metadata.get_tier_settings`):
        ``max_pes_per_hit_per_sipm`` (``max_pes_per_hit_combined`` for the
        combined map), ``store_expected_pes``,
        ``photoelectron_resolution_sigma`` and ``time_resolution_in_ns``.
    rows_per_chunk
        Number of scintillator hits to read at once. A window is never split
        between chunks, so a chunk can hold more.
    perf_block
        Context manager of a :func:`reboost.profile.make_profiler`, used to time the
        reboost calls.

    Returns
    -------
    dict[str, NumdetStats]
        For each light source with an optical map and hits in `windows`, the
        energy depositions outside the map or in bins without detection
        probability. Empty, with nothing written, if `n_windows` is zero.
    """
    if n_windows == 0:
        return {}

    max_pes_per_hit = (
        settings["max_pes_per_hit_combined"]
        if sipm == "all"
        else settings["max_pes_per_hit_per_sipm"]
    )
    store_expected_pes = settings["store_expected_pes"]

    with perf_block("load_optmap()"):
        optmaps = {
            name: reboost.spms.load_optmap(m, sipm) if isinstance(m, (str, Path)) else m
            for name, m in optmaps.items()
            if name in windows
        }

    # constant for this SiPM, so resolve them once instead of per chunk
    map_scaling = {}
    for name, source in sources.items():
        map_scaling[name] = resolve_map_scaling(source["optmap_scaling_factor"], sipm)
        msg = f"using optical map scaling factor {map_scaling[name]} for {name}/{sipm}"
        log.debug(msg)

    stats = {
        name: NumdetStats()
        for name, m in optmaps.items()
        if isinstance(m, OptmapForConvolve)
    }

    # hits of all sources per window, as cumulative counts
    cum_rows = np.cumsum(
        np.bincount(
            np.concatenate([w for _, w in windows.values()]), minlength=n_windows
        )
    )

    w_start = 0
    while w_start < n_windows:
        base = cum_rows[w_start - 1] if w_start > 0 else 0
        w_end = max(
            w_start + 1,
            int(np.searchsorted(cum_rows, base + rows_per_chunk, side="right")),
        )
        w_end = min(w_end, n_windows)

        pe_times, expected_pes, t0, evtid, window_idx = ([] for _ in range(5))
        for name, (first_row, w) in windows.items():
            i, j = np.searchsorted(w, [w_start, w_end])
            if i == j:
                continue

            stp_chunk = lh5.read(
                f"stp/{sources[name]['table']}",
                stp_file,
                start_row=first_row + i,
                n_rows=j - i,
            )
            # a hit at the cap has max_pes_per_hit photoelectrons, so its
            # window is flagged as saturated after pooling
            if name in stats:
                times, _, expected, chunk_stats = optmap_to_pes(
                    stp_chunk.view_as("ak"),
                    optmaps[name],
                    sipm,
                    material=sources[name]["material"],
                    map_scaling=map_scaling[name],
                    max_pes_per_hit=max_pes_per_hit,
                    store_expected_pes=store_expected_pes,
                    perf_block=perf_block,
                )
                stats[name] += chunk_stats
            else:
                times, expected = light_model_to_pes(
                    stp_chunk,
                    optmaps[name],
                    material=sources[name]["material"],
                    map_scaling=map_scaling[name],
                    max_pes_per_hit=max_pes_per_hit,
                    perf_block=perf_block,
                )

            pe_times.append(ak.without_parameters(times))
            if store_expected_pes:
                expected_pes.append(np.asarray(expected, dtype=np.float64))
            hits = reboost.init_hit_table(stp_chunk)
            t0.append(hits.t0.nda)
            evtid.append(hits.evtid.nda)
            window_idx.append(w[i:j] - w_start)

        with perf_block("pool_photoelectrons()"):
            window_table, pooled, is_saturated, expected = pool_photoelectrons(
                ak.concatenate(pe_times),
                np.concatenate(expected_pes) if store_expected_pes else None,
                np.concatenate(t0),
                np.concatenate(evtid),
                np.concatenate(window_idx),
                w_end - w_start,
                max_pes_per_hit,
            )

        _write_chunk(
            window_table,
            pooled,
            is_saturated,
            expected,
            out_file,
            sipm,
            sipm_uid,
            runid,
            usability,
            settings,
            perf_block,
        )
        w_start = w_end

    return stats
