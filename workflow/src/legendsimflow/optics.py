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

The opt tier holds, for each SiPM channel, the photoelectrons detected in each
hit of the liquid argon. Photoelectrons come either from the convolution of the
energy depositions with an optical map (``light_source: optmap``), or from the
optical photons tracked by Geant4 up to the SiPMs (``light_source:
tracked_photons``). Times are in nanoseconds, amplitudes in photoelectrons.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Iterable, Mapping
from contextlib import AbstractContextManager, nullcontext
from pathlib import Path

import awkward as ak
import numpy as np
import reboost
import reboost.spms
from lgdo import Array, Table, VectorOfVectors
from numpy.typing import NDArray
from reboost.optmap.convolve import NumdetStats, OptmapForConvolve

from legendsimflow import metadata as mutils

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
    map_scaling: float,
    max_pes_per_hit: int,
    store_expected_pes: bool,
    perf_block: Callable[[str], AbstractContextManager] = nullcontext,
) -> tuple[ak.Array, NDArray[np.bool_], NDArray[np.floating] | None, NumdetStats]:
    """Sample the photoelectrons of the energy depositions in the argon.

    The number of scintillation photons emitted by each step follows from its
    energy deposition and particle type. The optical map gives the probability
    that a photon emitted at the step position reaches `sipm`. The emission
    times follow the argon scintillation time profile.

    Parameters
    ----------
    chunk
        Rows of the argon step table (fields ``edep``, ``particle``, ``xloc``,
        ``yloc``, ``zloc`` and ``time``, one list per hit).
    optmap
        The optical map, already loaded with :func:`reboost.spms.load_optmap`.
    sipm
        Name of the SiPM, or ``"all"`` for the combined map.
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
            chunk.edep, chunk.particle, "lar"
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
            nr_pe, chunk.particle, chunk.time, "lar"
        )
        # the photoelectron_times() processor does not guarantee time ordering
        pe_times = ak.sort(pe_times, axis=-1)

    return pe_times, is_saturated, expected_pes, stats


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
    """Assemble the opt tier table of one SiPM for a chunk of hits.

    The table carries the event fields of `stp_chunk` (see
    :func:`reboost.init_hit_table`) and, per hit: ``dt`` (photoelectron times
    relative to the hit time ``t0``, in nanoseconds), ``energy`` (the
    amplitudes, in photoelectrons), ``is_saturated``, ``expected_pes`` (only
    if given), and the ``period``, ``run`` and encoded ``usability`` of the
    data run the hit belongs to.
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


def write_sipm_photoelectrons(
    chunks: Iterable[Table],
    out_file: str | Path,
    sipm: str,
    sipm_uid: int,
    runid: str,
    usability: str,
    *,
    settings: Mapping,
    optmap: str | Path | OptmapForConvolve | None = None,
    perf_block: Callable[[str], AbstractContextManager] = nullcontext,
) -> NumdetStats | None:
    """Write the photoelectrons of one SiPM to the opt tier file.

    Each row of `chunks` gives one row of the output table, appended to
    ``hit/<sipm>`` in `out_file` (``hit/spms`` for the combined map). With
    ``light_source: optmap``, `chunks` runs over the argon step table and
    `optmap` is the map, or its file name, which is then loaded here and
    released on return (the per-channel maps are too large to keep all of them
    in memory). With ``light_source: tracked_photons``, `chunks` runs over the
    step table of the SiPM.

    Parameters
    ----------
    chunks
        The step table rows to process, as returned by
        :class:`~lh5.io.iterator.LH5Iterator`.
    out_file
        The opt tier file.
    sipm
        Name of the SiPM, or ``"all"`` for the combined map.
    sipm_uid
        UID of the output table, the one of the SiPM or, for the combined map,
        of the argon volume.
    runid
        The data run the rows belong to, as in the simulation statistics
        partition file.
    usability
        Usability of the SiPM in `runid`.
    settings
        The ``opt`` tier settings (see :func:`.metadata.get_tier_settings`):
        ``light_source``, ``optmap_scaling_factor``,
        ``max_pes_per_hit_per_sipm`` (``max_pes_per_hit_combined`` for the
        combined map), ``store_expected_pes``,
        ``photoelectron_resolution_sigma`` and ``time_resolution_in_ns``.
    optmap
        Optical map or its file name. Required with ``light_source: optmap``.
    perf_block
        Context manager of a :func:`reboost.profile.make_profiler`, used to time the
        reboost calls.

    Returns
    -------
    NumdetStats or None
        With an optical map, the energy depositions outside the map or in
        bins without detection probability, summed over `chunks`. ``None``
        with tracked photons.
    """
    light_source = settings["light_source"]
    max_pes_per_hit = (
        settings["max_pes_per_hit_combined"]
        if sipm == "all"
        else settings["max_pes_per_hit_per_sipm"]
    )

    if light_source == "optmap":
        with perf_block("load_optmap()"):
            if not isinstance(optmap, OptmapForConvolve):
                optmap = reboost.spms.load_optmap(optmap, sipm)

        # constant for this SiPM, so resolve it once instead of per chunk
        map_scaling = resolve_map_scaling(settings["optmap_scaling_factor"], sipm)
        msg = f"using optical map scaling factor {map_scaling} for {sipm}"
        log.debug(msg)

        stats = NumdetStats()
    else:
        stats = None

    for stp_chunk in chunks:
        chunk = stp_chunk.view_as("ak")

        if light_source == "tracked_photons":
            pe_times, is_saturated = tracked_photons_to_pes(chunk, max_pes_per_hit)
            expected_pes = None
        else:
            pe_times, is_saturated, expected_pes, chunk_stats = optmap_to_pes(
                chunk,
                optmap,
                sipm,
                map_scaling=map_scaling,
                max_pes_per_hit=max_pes_per_hit,
                store_expected_pes=settings["store_expected_pes"],
                perf_block=perf_block,
            )
            stats += chunk_stats

        pe_times, pe_amps = apply_photoelectron_resolution(
            pe_times,
            settings["photoelectron_resolution_sigma"],
            settings["time_resolution_in_ns"],
            perf_block=perf_block,
        )

        with perf_block("write_hit_table_chunk()"):
            out_table = build_opt_table(
                stp_chunk,
                pe_times,
                pe_amps,
                is_saturated,
                expected_pes,
                runid,
                usability,
            )
            reboost.write_hit_table_chunk(
                out_table,
                "hit/" + ("spms" if sipm == "all" else sipm),
                out_file,
                uid=sipm_uid,
            )

    return stats
