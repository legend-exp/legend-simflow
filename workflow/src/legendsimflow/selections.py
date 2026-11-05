# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>
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
"""Event selections of the pdf tier, for simulations and LEGEND-200 data.

The pdf tier histograms the multiplicity-1 events in categories, one level
per cut: the LAr veto, the A/E and, in simulations, the energy in the FCCD
layer. Each event is in exactly one category of each level, so any selection
is a sum of categories (see :func:`.pdf.read_pdf`). The data is selected with
the same functions, so that data and PDFs compare bin by bin. See
:ref:`selections` in the user manual for how to use them.

The simulations (`cvt` tier) and the data (`evt` or `pet` tier) store the
same information, with a few different names that the functions resolve. The
selections are:

- data cleaning: no forced trigger, pulser or muon event, HPGe hits of good
  quality, see :func:`is_clean_event`. Data only;
- hit: HPGe energy above the `evt` tier threshold ``geds_energy_thr_kev``,
  25 keV in the data production and by default here, see
  :func:`is_hpge_hit_above_thr`. The functions below find the hits themselves;
- multiplicity: number of hits in the event, all in ON detectors
  (``geds.quality.is_good_channel``), see :func:`has_hpge_multiplicity`. The
  functions of the A/E and FCCD cuts read the one hit of a multiplicity-1
  event, as :func:`one_hit` does;
- LAr veto: the SiPM multiplicity and photoelectron sum against the
  ``lar_veto_multiplicity_thr`` and ``lar_veto_energy_sum_pe_thr`` settings of
  the `evt` tier, see :func:`is_spms_coincident` and :func:`lar_categories`;
- A/E: a valid A/E needs a usable PSD (``psd.is_good``) and a modeled A/E:
  simulated (``has_aoe``) or, in data, a detector whose A/E the simulation
  models (``is_modeled``, see :func:`attach_modelable_flag`). Single-site is
  above the low-side cut, high A/E above the high-side cut (pulse-shape
  library only). See :func:`aoe_flags` and :func:`aoe_categories`;
- FCCD: energy deposited in the FCCD layer, from the surface steps of the
  simulation, see :func:`fccd_categories`.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping

import awkward as ak
import numpy as np
from dbetto.utils import load_dict

from . import SimflowConfig, nersc, patterns
from .exceptions import SimflowConfigError
from .metadata import get_sanitized_fccd, parse_runid, runinfo

log = logging.getLogger(__name__)


def is_clean_event(evt: ak.Array) -> ak.Array:
    """Return the mask of the data events to analyse.

    Not a forced trigger, not in coincidence with the pulser or the muon veto,
    and with HPGe hits of good quality (``geds.quality.is_bb_like``). The
    simulation has none of these, so it needs no cleaning.

    Parameters
    ----------
    evt
        events of the `evt` or `pet` tier of the data.
    """
    return (
        ~evt.trigger.is_forced
        & ~evt.coincident.puls
        & ~evt.coincident.muon
        & ~evt.coincident.muon_offline
        & evt.geds.quality.is_bb_like
    )


def is_hpge_hit_above_thr(evt: ak.Array, energy_thr_kev: float = 25) -> ak.Array:
    """Return the mask of the HPGe hits with ``geds.energy`` above `energy_thr_kev` (keV).

    The threshold is applied again to the stored energies, so that the
    multiplicity matches the stored float32 energies and a threshold above the
    one of the `evt` tier can be chosen. The data production stores only the
    hits above 25 keV, so a lower threshold cannot be applied to the data.

    Parameters
    ----------
    evt
        events of the `evt` tier.
    energy_thr_kev
        energy threshold in keV. The default is the threshold of the data
        production and of the `evt` tier (``geds_energy_thr_kev``).
    """
    return evt.geds.energy > energy_thr_kev


def has_hpge_multiplicity(
    evt: ak.Array,
    min_mult: int,
    max_mult: int | None = None,
    energy_thr_kev: float = 25,
) -> ak.Array:
    """Return the mask of the events with `min_mult` to `max_mult` hits, all in ON detectors.

    `max_mult` defaults to `min_mult`. An event with an AC detector
    (``geds.quality.is_good_channel`` false) above threshold is never selected.

    Parameters
    ----------
    evt
        events of the `evt` tier.
    min_mult
        smallest number of hits, included.
    max_mult
        largest number of hits, included.
    energy_thr_kev
        hit energy threshold in keV, see :func:`is_hpge_hit_above_thr`.
    """
    if max_mult is None:
        max_mult = min_mult
    is_hit = is_hpge_hit_above_thr(evt, energy_thr_kev)
    mult = ak.sum(is_hit, axis=-1)
    all_on = ak.all(evt.geds.quality.is_good_channel[is_hit], axis=-1)
    return (mult >= min_mult) & (mult <= max_mult) & all_on


def one_hit(evt: ak.Array, energy_thr_kev: float = 25) -> ak.Array:
    """Return the HPGe record of the one hit of each event: the per-hit ``geds`` fields.

    For multiplicity-1 events, see :func:`has_hpge_multiplicity`: the one hit
    above threshold. ``one_hit(evt).energy`` is the energy of the event. The
    per-event fields of ``geds``, like ``multiplicity``, are left out.

    Parameters
    ----------
    evt
        multiplicity-1 events of the `evt` tier.
    energy_thr_kev
        hit energy threshold in keV, see :func:`is_hpge_hit_above_thr`.
    """
    return _one_hit(evt.geds, is_hpge_hit_above_thr(evt, energy_thr_kev))


def _one_hit(table: ak.Array, is_hit: ak.Array) -> ak.Array | None:
    """Return the record of the per-hit fields of `table` at the one selected hit."""
    fields = {}
    for name in table.fields:
        field = table[name]
        if field.fields:
            hit = _one_hit(field, is_hit)
            if hit is not None:
                fields[name] = hit
        # a field left out of the LH5 field mask is read as empty lists
        elif field.ndim > 1 and ak.all(ak.num(field) == ak.num(is_hit)):
            fields[name] = field[is_hit][:, 0]
    return ak.zip(fields, depth_limit=1) if fields else None


def aoe_model(hit_settings: Mapping) -> str | None:
    """Return the A/E model simulated with the `hit` tier settings.

    ``pulse_lib`` with ``simulate_psd_with_psl``, else ``single_temp`` with
    ``simulate_psd``, else `None`.
    """
    if hit_settings["simulate_psd_with_psl"]:
        return "pulse_lib"
    if hit_settings["simulate_psd"]:
        return "single_temp"
    return None


def category_levels(
    has_spms_coinc: bool, aoe_model: str | None, fccd_thr: float | None
) -> list[dict[str, str]]:
    """Return the event categories (name -> description), one dict per level.

    The names are the keys of the event masks of :func:`lar_categories`,
    :func:`aoe_categories` and :func:`fccd_categories`. A level is left out when
    its data is missing.

    Parameters
    ----------
    has_spms_coinc
        whether the SiPM data of the LAr veto is available.
    aoe_model
        ``pulse_lib``, ``single_temp`` or `None` without simulated PSD.
    fccd_thr
        energy threshold in keV of the FCCD tag, or `None` without it.
    """
    levels = []
    if has_spms_coinc:
        levels.append(
            {
                "is_spms_coincident": "light seen by the SiPMs (fails the LAr veto)",
                "not_spms_coincident": "no light seen by the SiPMs (passes the LAr veto)",
            }
        )
    if aoe_model is not None:
        aoe_level = {}
        # the single template has no high-side cut
        if aoe_model == "pulse_lib":
            aoe_level["is_high_aoe"] = "A/E above the high-side cut"
        aoe_level["is_ss_aoe"] = "single-site A/E (passes the A/E cuts)"
        aoe_level["is_low_aoe"] = "A/E below the low-side cut (multi-site)"
        aoe_level["is_invalid_aoe"] = "no valid or simulated A/E"
        levels.append(aoe_level)
    if fccd_thr is not None:
        levels.append(
            {
                "is_fccd_tagged": f"more than {fccd_thr} keV deposited in the FCCD layer",
                "not_fccd_tagged": f"at most {fccd_thr} keV deposited in the FCCD layer",
            }
        )
    return levels


def is_spms_coincident(
    evt: ak.Array, multiplicity_thr: int = 4, energy_sum_thr_pe: float = 4
) -> np.ndarray:
    """Return the LAr veto flag of the events: light seen by the SiPMs.

    The event fails the LAr veto when at least `multiplicity_thr` SiPM channels
    see a pulse (``spms.multiplicity``) or when the pulses sum to
    `energy_sum_thr_pe` photoelectrons or more (``spms.energy_sum``), the
    ``lar_veto_multiplicity_thr`` and ``lar_veto_energy_sum_pe_thr`` settings of
    the `evt` tier. This is ``coincident.spms`` of the simulated `evt` tier. In
    data the flag is computed here, not read from the data production, which
    counts the physical pulses above its noise threshold in the trigger
    coincidence window.

    Parameters
    ----------
    evt
        events of the `evt` tier.
    multiplicity_thr
        number of channels with a pulse that fails the veto. The default is
        the one of the LEGEND-200 analysis.
    energy_sum_thr_pe
        summed photoelectrons that fail the veto. The default is the one of
        the LEGEND-200 analysis.
    """
    spms = evt.spms
    return ak.to_numpy(
        (spms.multiplicity >= multiplicity_thr) | (spms.energy_sum >= energy_sum_thr_pe)
    )


def lar_categories(is_spms_coincident: np.ndarray) -> dict[str, np.ndarray]:
    """Return the event masks of the LAr veto categories.

    The categories are defined by the LAr veto flag of
    :func:`is_spms_coincident`, `True` when the SiPMs see enough light:

    - ``is_spms_coincident``: light seen, the event fails the LAr veto;
    - ``not_spms_coincident``: no light seen, the event passes the LAr veto.

    Parameters
    ----------
    is_spms_coincident
        the LAr veto flag of each event.
    """
    spms = np.asarray(is_spms_coincident, dtype=bool)
    return {"is_spms_coincident": spms, "not_spms_coincident": ~spms}


def attach_modelable_flag(
    evt: ak.Array, config: SimflowConfig, is_modelable: Mapping | None = None
) -> ak.Array:
    """Return the data events with ``geds.psd.is_modeled``: the hit is in a detector whose A/E the simulation models.

    The detectors are those flagged in the modelable-HPGe cache of the `par`
    step (``pars/detinfo/is_modelable.yaml`` of the production) for the run of
    each event, ``trigger.period`` and ``trigger.run``. Their rawids come from
    the channel map of the Simflow metadata at the run start, the one the
    simulated geometry was built from. A run of the events missing from the
    cache is an error.

    Parameters
    ----------
    evt
        events of the `evt` or `pet` tier of the data.
    config
        Simflow configuration, see :func:`.utils.load_simflow_config`.
    is_modelable
        the cache (run id -> detector name -> flag) to use instead of the one
        of the production.
    """
    if is_modelable is None:
        is_modelable = load_dict(
            nersc.dvs_ro(config, patterns.detinfo_filename(config, "is_modelable"))
        )
    runids = {}
    for runid in is_modelable:
        _, period, run, _ = parse_runid(runid)
        runids.setdefault((period, run), []).append(runid)

    rawid = evt.geds.rawid
    period, run = (
        ak.to_numpy(ak.flatten(ak.broadcast_arrays(evt.trigger[f], rawid)[0]))
        for f in ("period", "run")
    )
    flat_rawid = ak.to_numpy(ak.flatten(rawid))
    is_modeled = np.zeros(len(flat_rawid), dtype=bool)
    for p, r in np.unique(np.stack([period, run], axis=1), axis=0):
        if len(runids.get((p, r), [])) != 1:
            msg = (
                f"period {p} run {r} of the data matches "
                f"{runids.get((p, r), [])} in the modelable-HPGe cache, need one run"
            )
            raise SimflowConfigError(msg)
        runid = runids[p, r][0]
        meta = config.metadata
        chmap = meta.channelmap(runinfo(meta, runid).start_key)
        rawids = [int(chmap[n].daq.rawid) for n, v in is_modelable[runid].items() if v]
        is_modeled |= (period == p) & (run == r) & np.isin(flat_rawid, rawids)
    is_modeled = ak.unflatten(is_modeled, ak.num(rawid))
    return ak.with_field(evt, is_modeled, ["geds", "psd", "is_modeled"])


def aoe_flags(
    evt: ak.Array, aoe_model: str, energy_thr_kev: float = 25
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Return the A/E flags of the one hit of each event: valid, single-site, high.

    The flags come from ``geds.psd``. The data (`evt` or `pet` tier) is
    recognized by its ``psd.low_aoe`` table:

    - valid: ``psd.is_good`` (valid PSD in data) and a modeled A/E:
      ``psd.<aoe_model>.has_aoe`` in the simulation, ``psd.is_modeled`` of
      :func:`attach_modelable_flag` in data;
    - single-site: the classifier is above the low-side cut,
      ``psd.<aoe_model>.is_single_site`` or ``psd.low_aoe.is_single_site`` in
      data;
    - high: the classifier is above the high-side cut,
      ``psd.pulse_lib.is_high_aoe`` or not ``psd.high_aoe.is_bulk`` in data.
      `None` with the single template, which has no high-side cut.

    Parameters
    ----------
    evt
        multiplicity-1 events of the `evt` tier, see :func:`one_hit`.
    aoe_model
        ``pulse_lib`` (pulse-shape library) or ``single_temp`` (single
        template), the model simulated in the `hit` tier.
    energy_thr_kev
        hit energy threshold in keV, see :func:`is_hpge_hit_above_thr`.
    """
    psd = one_hit(evt, energy_thr_kev).psd
    if "low_aoe" in psd.fields:
        if "is_modeled" not in psd.fields:
            msg = (
                "no geds.psd.is_modeled in the data, add it with attach_modelable_flag"
            )
            raise ValueError(msg)
        modeled = psd.is_modeled
        single_site = psd.low_aoe.is_single_site
        high = ~psd.high_aoe.is_bulk if aoe_model == "pulse_lib" else None
    else:
        aoe = psd[aoe_model]
        modeled = aoe.has_aoe
        single_site = aoe.is_single_site
        high = aoe.is_high_aoe if aoe_model == "pulse_lib" else None
    valid = ak.to_numpy(psd.is_good & modeled)
    single_site = ak.to_numpy(single_site)
    if high is not None:
        high = ak.to_numpy(high)
    return valid, single_site, high


def aoe_categories(
    valid: np.ndarray, single_site: np.ndarray, high: np.ndarray | None = None
) -> dict[str, np.ndarray]:
    """Return the event masks of the A/E categories, from the flags of one hit.

    The flags are those of :func:`aoe_flags`. In order of precedence:

    - ``is_invalid_aoe``: no valid A/E;
    - ``is_high_aoe``: above the high-side cut. Only with `high`, the single
      template has no high-side cut;
    - ``is_ss_aoe``: single-site, above the low-side cut;
    - ``is_low_aoe``: multi-site, below the low-side cut.

    Each event is in exactly one category.

    Parameters
    ----------
    valid
        whether the A/E of the event is valid.
    single_site
        whether the A/E is above the low-side cut.
    high
        whether the A/E is above the high-side cut, or `None` without
        high-side cut.
    """
    categories = {}
    if high is None:
        high = np.zeros(len(valid), dtype=bool)
    else:
        high = valid & high
        categories["is_high_aoe"] = high
    categories["is_ss_aoe"] = valid & ~high & single_site
    categories["is_low_aoe"] = valid & ~high & ~single_site
    categories["is_invalid_aoe"] = ~valid
    return categories


def pick_detector(evt: ak.Array, rawid: int, table: str = "geds") -> ak.Array:
    """Return the events with a hit in one detector, with the fields of that hit only.

    The per-hit fields of `table` (for example ``geds.energy`` and the steps in
    ``geds.surface``) lose the hit level: ``energy`` gets one value per event,
    ``surface.edep`` one list of steps per event. The fields with one value
    per event (for example ``geds.energy_sum``) and the other tables of `evt`
    are kept as they are, for the selected events. The hits are taken as
    stored, with no energy threshold.

    Parameters
    ----------
    evt
        events of the `evt` tier.
    rawid
        the detector, as in the ``rawid`` field of `table`.
    table
        the detector table, ``geds`` or ``spms``.

    Examples
    --------
    >>> v = pick_detector(evt, 1104000)
    >>> v.geds.energy  # one energy per event, in keV
    >>> v.geds.surface.distance_to_nplus  # one list of steps per event, in mm
    >>> v.coincident.spms  # one value per selected event
    """
    is_det = evt[table].rawid == rawid
    if ak.any(ak.sum(is_det, axis=-1) > 1):
        msg = f"detector {rawid} has more than one hit in an event"
        raise ValueError(msg)

    in_det = ak.any(is_det, axis=-1)
    evt = evt[in_det]
    is_det = is_det[in_det]

    def _pick(field):
        if field.fields:
            return ak.zip({k: _pick(field[k]) for k in field.fields}, depth_limit=1)
        if field.ndim == 1:
            return field
        return ak.flatten(field[is_det], axis=1)

    # packed: repeated computations on the steps (e.g. a PDF grid scan) run faster
    return ak.with_field(evt, ak.to_packed(_pick(evt[table])), table)


def fccd_categories(
    evt: ak.Array,
    fccd_mm: np.ndarray,
    fccd_thr_kev: float,
    energy_thr_kev: float = 25,
) -> dict[str, np.ndarray]:
    """Return the event masks of the FCCD categories, from the steps of one hit.

    The energy in the FCCD layer is the sum of the energy depositions
    ``geds.surface.edep`` (keV) with a distance to the n+ surface
    ``geds.surface.distance_to_nplus`` (mm) smaller than the FCCD of the
    detector. These are the energies as deposited, before the dead-layer
    model:

    - ``is_fccd_tagged``: more than `fccd_thr_kev` in the FCCD layer;
    - ``not_fccd_tagged``: at most `fccd_thr_kev` in the FCCD layer.

    Parameters
    ----------
    evt
        multiplicity-1 events of the `cvt` tier, see :func:`one_hit`.
    fccd_mm
        FCCD in mm of the detector of the hit, for each event.
    fccd_thr_kev
        energy threshold in keV of the tag.
    energy_thr_kev
        hit energy threshold in keV, see :func:`is_hpge_hit_above_thr`.
    """
    surface = one_hit(evt, energy_thr_kev).surface
    in_fccd = surface.distance_to_nplus < fccd_mm
    edep_fccd = ak.sum(surface.edep[in_fccd], axis=-1)
    tagged = ak.to_numpy(edep_fccd > fccd_thr_kev)
    return {"is_fccd_tagged": tagged, "not_fccd_tagged": ~tagged}


def fccd_by_uid(
    metadata: Mapping,
    name_to_uid: Mapping[str, int],
    surface_depth: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the sorted detector UIDs and their FCCD in mm.

    The FCCD is read from the detector metadata with
    :func:`.metadata.get_sanitized_fccd`. It is limited to `surface_depth`,
    the depth of the steps stored in ``geds.surface``, with a warning. The
    FCCD of the detector of each hit is then
    ``fccd_mm[np.searchsorted(uids, uid)]``.

    Parameters
    ----------
    metadata
        LEGEND metadata database.
    name_to_uid
        map of the detector names to their UIDs, e.g. ``detector_uids`` of
        the `cvt` tier.
    surface_depth
        ``surface_depth_in_mm`` of the `hit` tier settings.
    """
    names = sorted(name_to_uid, key=name_to_uid.get)
    uids = np.array([name_to_uid[n] for n in names], dtype=np.int64)
    fccd_mm = np.zeros(len(names))
    for i, name in enumerate(names):
        fccd_mm[i] = get_sanitized_fccd(metadata, name)
        # the steps deeper than surface_depth are not stored
        if fccd_mm[i] > surface_depth:
            log.warning(
                "FCCD of %s (%s mm) exceeds surface_depth_in_mm, using %s mm",
                name,
                fccd_mm[i],
                surface_depth,
            )
            fccd_mm[i] = surface_depth
    return uids, fccd_mm
