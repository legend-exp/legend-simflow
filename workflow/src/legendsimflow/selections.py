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
the same functions, so that data and PDFs compare bin by bin.

The simulations (`cvt` tier) and the data (`evt` or `pet` tier) store the
same information under different names. The selections are:

- hit: HPGe energy above the `evt` tier threshold ``geds_energy_thr_kev``, in
  an ON detector (``geds.quality.is_good_channel``), see :func:`select_hits`;
- multiplicity: number of hits in the event;
- data cleaning: no forced trigger, pulser or muon event, HPGe hits of good
  quality, see :func:`is_clean_data_event`. Data only;
- LAr veto: ``coincident.spms``, see :func:`lar_categories`;
- A/E: a valid A/E needs a usable PSD (``psd.is_good``) and a modeled A/E:
  simulated (``has_aoe``) or, in data, a detector whose A/E the simulation
  models. Single-site is above the low-side cut, high A/E above the
  high-side cut (pulse-shape library only). See :func:`sim_aoe_flags`,
  :func:`data_aoe_flags` and :func:`aoe_categories`;
- FCCD: energy deposited in the FCCD layer, from the surface steps of the
  simulation, see :func:`fccd_categories`.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence

import awkward as ak
import numpy as np

from .metadata import get_sanitized_fccd

log = logging.getLogger(__name__)


def select_hits(
    evt: ak.Array, energy_thr_kev: float
) -> tuple[ak.Array, ak.Array, ak.Array]:
    """Return the hits above threshold, the multiplicity and whether all hits are in ON detectors.

    The multiplicity is counted here, not read from ``geds.multiplicity``, so
    that it matches the stored float32 energies.

    Parameters
    ----------
    evt
        events of the `evt` tier.
    energy_thr_kev
        energy threshold in keV.
    """
    is_hit = evt.geds.energy > energy_thr_kev
    multiplicity = ak.sum(is_hit, axis=-1)
    is_good = ak.all(evt.geds.quality.is_good_channel[is_hit], axis=-1)
    return is_hit, multiplicity, is_good


def one_hit(field: ak.Array, is_hit: ak.Array) -> ak.Array:
    """Return the field of the one selected hit of each event."""
    return field[is_hit][:, 0]


def is_clean_data_event(evt: ak.Array) -> ak.Array:
    """Return the mask of the data events to analyse.

    Not a forced trigger, not in coincidence with the pulser or the muon veto,
    and with HPGe hits of good quality (``geds.quality.is_bb_like``).

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
        whether the LAr veto flag ``coincident.spms`` is available.
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


def lar_categories(evt: ak.Array) -> dict[str, np.ndarray]:
    """Return the event masks of the LAr veto categories.

    The categories are defined by ``coincident.spms``, the LAr veto flag of
    the `evt` tier, which is `True` when the SiPMs see enough light:

    - ``is_spms_coincident``: light seen, the event fails the LAr veto;
    - ``not_spms_coincident``: no light seen, the event passes the LAr veto.

    Parameters
    ----------
    evt
        events of the `evt` tier.
    """
    spms = ak.to_numpy(evt.coincident.spms)
    return {"is_spms_coincident": spms, "not_spms_coincident": ~spms}


def sim_aoe_flags(
    evt: ak.Array, is_hit: ak.Array, aoe_model: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Return the A/E flags of the one hit of each simulated event.

    The flags come from ``geds.psd`` and its ``<aoe_model>`` subtable:

    - valid: ``psd.is_good`` (valid PSD in data) and ``has_aoe`` (A/E
      simulated);
    - single-site: ``is_single_site``, the classifier is above the low-side
      cut;
    - high: ``is_high_aoe``, the classifier is above the high-side cut. `None`
      with the single template, which has no high-side cut.

    Parameters
    ----------
    evt
        events of the `cvt` tier, with one hit selected by `is_hit` in each.
    is_hit
        selection of the hit of each event.
    aoe_model
        ``pulse_lib`` (pulse-shape library) or ``single_temp`` (single
        template).
    """
    psd = evt.geds.psd
    aoe = psd[aoe_model]
    valid = ak.to_numpy(one_hit(psd.is_good & aoe.has_aoe, is_hit))
    single_site = ak.to_numpy(one_hit(aoe.is_single_site, is_hit))
    high = None
    if aoe_model == "pulse_lib":
        high = ak.to_numpy(one_hit(aoe.is_high_aoe, is_hit))
    return valid, single_site, high


def data_aoe_flags(
    evt: ak.Array, is_hit: ak.Array, aoe_model: str, modeled_rawids: Sequence[int]
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Return the A/E flags of the one hit of each data event, as :func:`sim_aoe_flags`.

    The flags come from ``geds.psd`` of the `evt` or `pet` tier:

    - valid: ``psd.is_good`` and the hit is in a detector of `modeled_rawids`;
    - single-site: ``low_aoe.is_single_site``;
    - high: not ``high_aoe.is_bulk``. `None` with the single template.

    Parameters
    ----------
    evt
        events of the `evt` or `pet` tier, with one hit selected by `is_hit`
        in each.
    is_hit
        selection of the hit of each event.
    aoe_model
        ``pulse_lib`` or ``single_temp``, the model simulated in the `hit`
        tier.
    modeled_rawids
        ``rawid`` of the detectors whose A/E the simulation models.
    """
    psd = evt.geds.psd
    rawid = ak.to_numpy(one_hit(evt.geds.rawid, is_hit))
    valid = ak.to_numpy(one_hit(psd.is_good, is_hit)) & np.isin(rawid, modeled_rawids)
    single_site = ak.to_numpy(one_hit(psd.low_aoe.is_single_site, is_hit))
    high = None
    if aoe_model == "pulse_lib":
        high = ~ak.to_numpy(one_hit(psd.high_aoe.is_bulk, is_hit))
    return valid, single_site, high


def aoe_categories(
    valid: np.ndarray, single_site: np.ndarray, high: np.ndarray | None = None
) -> dict[str, np.ndarray]:
    """Return the event masks of the A/E categories, from the flags of one hit.

    The flags are those of :func:`sim_aoe_flags` and :func:`data_aoe_flags`.
    In order of precedence:

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


def fccd_categories(
    evt: ak.Array, is_hit: ak.Array, fccd_mm: np.ndarray, energy_thr_kev: float
) -> dict[str, np.ndarray]:
    """Return the event masks of the FCCD categories, from the steps of one hit.

    The energy in the FCCD layer is the sum of the energy depositions
    ``geds.surface.edep`` (keV) with a distance to the n+ surface
    ``geds.surface.distance_to_nplus`` (mm) smaller than the FCCD of the
    detector. These are the energies as deposited, before the dead-layer
    model:

    - ``is_fccd_tagged``: more than `energy_thr_kev` in the FCCD layer;
    - ``not_fccd_tagged``: at most `energy_thr_kev` in the FCCD layer.

    Parameters
    ----------
    evt
        events of the `cvt` tier, with one hit selected by `is_hit` in each.
    is_hit
        selection of the hit of each event.
    fccd_mm
        FCCD in mm of the detector of the hit, for each event.
    energy_thr_kev
        energy threshold in keV.
    """
    surface = evt.geds.surface
    in_fccd = one_hit(surface.distance_to_nplus, is_hit) < fccd_mm
    edep_fccd = ak.sum(one_hit(surface.edep, is_hit)[in_fccd], axis=-1)
    tagged = ak.to_numpy(edep_fccd > energy_thr_kev)
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
