"""
Module to implement the dataset preparation and average ("superpulse") construction to characterize the pulse shape response of HPGe detectors.

This is an important step in tuning the pulse shape discrimination (PSD) simulations

Conventions
-----------
- All event-level data is handled as ``ak.Array`` (awkward-array).
- Waveform data on disk is read via ``lgdo.lh5.read()`` and converted to
  ``ak.Array`` for in-memory processing.
- All times are in nanoseconds unless otherwise stated.
- All energies are in keV unless otherwise stated.
- The ``tab_map`` convention follows the format:
  ``{detector_name: rawid}`` e.g. ``{"V03422A": 1084803}``

Pipeline order (per file, inside the orchestrator loop)
--------------------------------------------------------
Steps 1-7 are called once per file inside ``build_superpulse_for_slice``,
which accumulates results across files until enough waveforms are available.

1. ``read_evt_data``               - read EVT tier for one file (no DSP fields)
2. ``perform_data_selection``      - quality + PSD cuts on EVT fields only
3. ``select_detector_events``      - filter to one detector
4. ``add_dsp_pars_to_evt``         - attach DSP fields for one detector only
5. ``select_data_in_slice``        - filter to one energy-drift-time slice
6. ``get_charge_and_current_wfs_for_slice`` - run production DSP chain via
                                      WaveformBrowser, return aligned charge
                                      and current waveforms for all slice events

Steps 7-11 run once, after enough waveforms have been accumulated:

7.  ``compute_superpulse``          - preliminary average (charge + current)
8.  ``compute_chi2_vs_superpulse``  - reduced chi2 of each wf vs superpulse
9.  ``apply_chi2_cut``              - retain golden waveforms only
10. ``compute_superpulse``          - final average on golden waveforms

If after step 9 the number of golden waveforms is below ``n_target_wfs``,
the orchestrator returns to step 1 and continues reading files from where
it left off, re-running steps 7-10 on the enlarged accumulated set, until
either the target is reached or all files are exhausted.

Step 11 runs once per detector, after all slices are done:

11. ``write_superpulses_to_lh5``   - write all slices for one detector to LH5
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import awkward as ak
import numpy as np
from lgdo import Array as LGDOArray
from lgdo import Scalar as LGDOScalar
from lgdo import Struct as LGDOStruct
from lgdo import lh5

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Slice:
    """
    Defines a 2D slice in the energy-drift-time space.

    Parameters
    ----------
    energy_range : tuple[float, float]
        Lower and upper bounds of the energy slice, in keV.
        Example: ``(1500.0, 2000.0)``
    drift_time_range : tuple[float, float]
        Lower and upper bounds of the drift time slice, in ns.
        Example: ``(900.0, 1100.0)``

    """

    energy_range: tuple[float, float]
    drift_time_range: tuple[float, float]

    @property
    def drift_time_center(self) -> float:
        """Center of the drift time range, in ns."""
        return (self.drift_time_range[0] + self.drift_time_range[1]) / 2.0

    @property
    def energy_center(self) -> float:
        """Center of the energy range, in keV."""
        return (self.energy_range[0] + self.energy_range[1]) / 2.0

    def __str__(self) -> str:
        return (
            f"Slice(E=[{self.energy_range[0]}, {self.energy_range[1]}] keV, "
            f"dt=[{self.drift_time_range[0]}, {self.drift_time_range[1]}] ns)"
        )


class Superpulse:
    """
    Holds the average charge and current waveforms for one slice.

    One ``Superpulse`` instance is produced per slice per detector, after the
    full preprocessing and chi2 self-similarity cut. It carries both waveforms
    together with the metadata needed to identify the slice, write to LH5, and
    perform the subsequent electronics parameter optimisation.

    Parameters
    ----------
    charge_wf : np.ndarray
        Average normalised charge waveform, shape ``(n_charge_samples,)``.
        Amplitude is dimensionless (ADC / cuspEmax, normalised to 1 at
        plateau).
    current_wf : np.ndarray
        Average current waveform after MWA, shape ``(n_current_samples,)``.
        Units: ``(ADC / cuspEmax) / ns * dt_data``, matching the convention
        in ``psl.py``.
    charge_time_axis : np.ndarray
        Time axis for the charge waveform in ns, shape ``(n_charge_samples,)``,
        aligned so that ``tp_aoe_max = 0``.
    current_time_axis : np.ndarray
        Time axis for the current waveform in ns, shape ``(n_current_samples,)``,
        aligned so that ``tp_aoe_max = 0``.
    slice : Slice
        The energy-drift-time slice this superpulse represents.
    detector : str
        Detector name, e.g. ``"V03422A"``.
    n_events_preliminary : int
        Number of waveforms used to build the preliminary superpulse (before
        the chi2 cut).
    n_events_final : int
        Number of waveforms surviving the chi2 cut, used to build this
        superpulse.

    Examples
    --------
    Build a superpulse and access its waveforms::

        sp = Superpulse(
            charge_wf=avg_charge,
            current_wf=avg_current,
            charge_time_axis=t_charge,
            current_time_axis=t_current,
            slice=Slice((1500., 2000.), (900., 1100.)),
            detector="V03422A",
            n_events_preliminary=120,
            n_events_final=98,
        )
        sp.charge_wf            # np.ndarray
        sp.current_wf           # np.ndarray
        sp.slice.drift_time_center  # 1000.0 ns

    Use as dict value (Slice is hashable)::

        superpulses: dict[Slice, Superpulse] = {}
        superpulses[sp.slice] = sp
    """

    def __init__(
        self,
        charge_wf: np.ndarray,
        current_wf: np.ndarray,
        charge_time_axis: np.ndarray,
        current_time_axis: np.ndarray,
        slice: Slice,
        detector: str,
        n_events_preliminary: int,
        n_events_final: int,
    ) -> None:
        if len(charge_wf) != len(charge_time_axis):
            msg = (
                f"charge_wf ({len(charge_wf)}) and charge_time_axis "
                f"({len(charge_time_axis)}) must have the same length."
            )
            raise ValueError(msg)
        if len(current_wf) != len(current_time_axis):
            msg = (
                f"current_wf ({len(current_wf)}) and current_time_axis "
                f"({len(current_time_axis)}) must have the same length."
            )
            raise ValueError(msg)
        if n_events_final > n_events_preliminary:
            msg = (
                f"n_events_final ({n_events_final}) cannot exceed "
                f"n_events_preliminary ({n_events_preliminary})."
            )
            raise ValueError(msg)

        self.charge_wf = charge_wf
        self.current_wf = current_wf
        self.charge_time_axis = charge_time_axis
        self.current_time_axis = current_time_axis
        self.slice = slice
        self.detector = detector
        self.n_events_preliminary = n_events_preliminary
        self.n_events_final = n_events_final

    def __repr__(self) -> str:
        return (
            f"Superpulse("
            f"detector={self.detector!r}, "
            f"slice={self.slice}, "
            f"n_events={self.n_events_final}/{self.n_events_preliminary}, "
            f"n_charge_samples={len(self.charge_time_axis)}, "
            f"n_current_samples={len(self.current_time_axis)})"
        )

    def to_lgdo(
        self,
    ):  ### Is this the best way to implement it? Can I use lgdo functions directly?
        """
        Serialise to an ``lgdo.Struct`` ready for writing to LH5.

        Returns
        -------
        lgdo.Struct
            Struct with the following fields:

            - ``charge_wf``            : ``lgdo.Array``, shape ``(n_charge_samples,)``
            - ``current_wf``           : ``lgdo.Array``, shape ``(n_current_samples,)``
            - ``charge_time_axis``     : ``lgdo.Array``, shape ``(n_charge_samples,)``,
              attrs ``{"units": "ns"}``
            - ``current_time_axis``    : ``lgdo.Array``, shape ``(n_current_samples,)``,
              attrs ``{"units": "ns"}``
            - ``dt_center``            : ``lgdo.Scalar``, drift time center [ns]
            - ``dt_lo``                : ``lgdo.Scalar``, drift time lower bound [ns]
            - ``dt_hi``                : ``lgdo.Scalar``, drift time upper bound [ns]
            - ``e_lo``                 : ``lgdo.Scalar``, energy lower bound [keV]
            - ``e_hi``                 : ``lgdo.Scalar``, energy upper bound [keV]
            - ``detector``             : ``lgdo.Scalar``, detector name string
            - ``n_events_preliminary`` : ``lgdo.Scalar``
            - ``n_events_final``       : ``lgdo.Scalar``
        """
        return LGDOStruct(
            {
                "charge_wf": LGDOArray(self.charge_wf),
                "current_wf": LGDOArray(self.current_wf),
                "charge_time_axis": LGDOArray(
                    self.charge_time_axis, attrs={"units": "ns"}
                ),
                "current_time_axis": LGDOArray(
                    self.current_time_axis, attrs={"units": "ns"}
                ),
                "dt_center": LGDOScalar(
                    self.slice.drift_time_center, attrs={"units": "ns"}
                ),
                "dt_lo": LGDOScalar(
                    self.slice.drift_time_range[0], attrs={"units": "ns"}
                ),
                "dt_hi": LGDOScalar(
                    self.slice.drift_time_range[1], attrs={"units": "ns"}
                ),
                "e_lo": LGDOScalar(self.slice.energy_range[0], attrs={"units": "keV"}),
                "e_hi": LGDOScalar(self.slice.energy_range[1], attrs={"units": "keV"}),
                "detector": LGDOScalar(self.detector),
                "n_events_preliminary": LGDOScalar(self.n_events_preliminary),
                "n_events_final": LGDOScalar(self.n_events_final),
            }
        )


# ---------------------------------------------------------------------------
# Dataset preparation (individual steps, one file at a time)
# ---------------------------------------------------------------------------


def read_evt_data(  ### Can / should this be merge with the data selection function? Code here is trivial
    evt_files: list[str],
) -> ak.Array:
    """
    Read EVT-tier LH5 files into a single concatenated awkward array.

    No DSP fields are attached here. This is intentionally a lightweight
    read so that ``perform_data_selection`` and ``select_detector_events``
    can reduce the array size before the DSP read in ``add_dsp_pars_to_evt``.

    In normal use this is called with a single-element list by the orchestrator
    ``build_superpulse_for_slice``, which manages the file loop externally.
    Passing a multi-file list is supported for testing or interactive use.

    Parameters
    ----------
    evt_files : list[str]
        Paths to one or more EVT-tier LH5 files. Results are concatenated
        along the event axis.

    Returns
    -------
    ak.Array
        Concatenated EVT data with shape ``(n_events,)``.
    """
    evt_arrays = []
    for evt_file in evt_files:
        evt_arrays.append(lh5.read_as("evt/", evt_file, library="ak"))

    if len(evt_arrays) == 1:
        return evt_arrays[0]
    return ak.concatenate(evt_arrays)


def perform_data_selection(
    evt_data: ak.Array,
    aoe_low_threshold: float = -3.0,
    aoe_high_threshold: float = 3.0,
) -> ak.Array:
    """
    Apply quality and PSD cuts to retain only valid single-site events.

    Two sequential sets of cuts are applied:

    1. Preliminary cuts: forced triggers, pulser events, muon coincidences
       (online and offline), bad channel quality, failing QC, multiplicity > 1,
       and insufficient LAr light (``spms.energy_sum <= 10``).

    2. PSD cuts: A/E lower and upper bounds on ``geds.psd.low_aoe.value``
       and ``geds.psd.high_aoe.value``. These cuts are needed to prevent
       the preliminary superpulse from being too contaminated for the subsequent
       chi2 self-similarity cut to be meaningful.

    Parameters
    ----------
    evt_data : ak.Array
        EVT data as returned by ``read_evt_data``.
    aoe_low_threshold : float, optional
        Lower bound on ``geds.psd.low_aoe.value``. Default: ``-3.0``.
    aoe_high_threshold : float, optional
        Upper bound on ``geds.psd.high_aoe.value``. Default: ``3.0``.

    Returns
    -------
    ak.Array
        Filtered event data with the same structure as the input.

    Notes
    -----
    The ``multiplicity == 1`` cut guarantees that all downstream ``geds.*``
    arrays have exactly one entry per event, simplifying all downstream
    indexing.
    """
    n_before = len(evt_data)

    # Preliminary cuts
    preliminary_mask = (
        ak.all(evt_data.geds.quality.is_good_channel, axis=-1)
        & ~evt_data.trigger.is_forced
        & ~evt_data.coincident.puls
        & ~evt_data.coincident.muon
        & ~evt_data.coincident.muon_offline
        & evt_data.geds.quality.is_bb_like
        & (evt_data.geds.multiplicity == 1)
        & (
            evt_data.spms.energy_sum > 10
        )  ### Is this actually needed? If not I could drop the spms entirely when reading evt files!
    )
    evt_data = evt_data[preliminary_mask]

    n_after_preliminary = len(evt_data)
    log.debug("preliminary cuts: %d -> %d events", n_before, n_after_preliminary)

    # PSD cuts
    psd_mask = ak.all(
        evt_data.geds.psd.low_aoe.value > aoe_low_threshold, axis=-1
    ) & ak.all(evt_data.geds.psd.high_aoe.value < aoe_high_threshold, axis=-1)
    evt_data = evt_data[psd_mask]

    log.debug("PSD cuts: %d -> %d events", n_after_preliminary, len(evt_data))
    return evt_data


def select_detector_events(
    evt_data: ak.Array,
    detector: str,
    tab_map: dict[str, int],
) -> ak.Array:
    """
    Filter event data to only events from a single detector.

    Called after ``perform_data_selection`` and before
    ``add_dsp_pars_to_evt``, so that the DSP read operates on the smallest
    possible array.

    Parameters
    ----------
    evt_data : ak.Array
        Event data after quality and PSD cuts, shape ``(n_events,)``.
    detector : str
        Name of the target detector, e.g. ``"V03422A"``.
    tab_map : dict[str, int]
        Mapping from detector name to rawid.

    Returns
    -------
    ak.Array
        Filtered event data containing only events where
        ``geds.rawid == tab_map[detector]``, shape ``(n_det_events,)``.

    Raises
    ------
    KeyError
        If ``detector`` is not found in ``tab_map``.
    """
    rawid = tab_map[detector]  # raises KeyError if missing
    mask = ak.any(evt_data.geds.rawid == rawid, axis=-1)
    det_evt_data = evt_data[mask]

    log.debug(
        "detector %s (rawid %d): %d -> %d events",
        detector,
        rawid,
        len(evt_data),
        len(det_evt_data),
    )
    return det_evt_data


def add_dsp_pars_to_evt(
    det_evt_data: ak.Array,
    dsp_file: str,
    tab_map: dict[str, int],
    fields: list[str],
) -> ak.Array:
    """
    Attach per-channel DSP parameters to a single-detector event array.

    Called after ``select_detector_events``, so ``det_evt_data`` already
    contains only events for one detector. The DSP read therefore involves
    only one channel: one ``lh5.read`` call per field.

    ``drift_time = tp_aoe_max - spms.first_t0`` is also computed and added
    as a top-level field, since it requires both a DSP field (``tp_aoe_max``)
    and an EVT field (``spms.first_t0``).

    Parameters
    ----------
    det_evt_data : ak.Array
        Single-detector event data as returned by ``select_detector_events``,
        shape ``(n_det_events,)``.
    dsp_file : str
        DSP-tier LH5 file path, corresponding to the same run segment as the
        EVT file used to produce ``det_evt_data``. Must be a single file
        because ``hit_idx`` values are row indices into this specific file.
    detector : str
        Name of the target detector, e.g. ``"V03422A"``.
    tab_map : dict[str, int]
        Mapping from detector name to rawid. Only the entry for ``detector``
        is used.
    fields : list[str]
        DSP fields to attach, e.g. ``["tp_aoe_max", "tp_0_est"]``.
        ``tp_aoe_max`` must be included for ``drift_time`` to be computed.

    Returns
    -------
    ak.Array
        ``det_evt_data`` with the requested DSP fields and ``drift_time``
        added as top-level fields. Each field has shape ``(1,)`` per event,
        consistent with the multiplicity == 1 guarantee.

    Raises
    ------
    ValueError
        If ``tp_aoe_max`` is not in ``fields`` (required for drift_time).
    """
    # Handle the empty-fields case: nothing to attach
    if not fields:
        log.debug("no DSP fields requested, skipping DSP read")
        return det_evt_data

    if "tp_aoe_max" not in fields:
        msg = (
            "'tp_aoe_max' must be included in fields "
            "(required for drift_time computation)"
        )
        raise ValueError(msg)

    # Read each DSP field and attach as a top-level field
    for field in fields:
        field_data = _read_dsp_field_for_channel(
            det_evt_data.geds.rawid,
            det_evt_data.geds.hit_idx,
            dsp_file,
            field,
            tab_map,
        )
        det_evt_data = ak.with_field(det_evt_data, field_data, field)

    # Compute drift time: tp_aoe_max - spms.first_t0
    det_evt_data = ak.with_field(
        det_evt_data,
        det_evt_data.tp_aoe_max - det_evt_data.spms.first_t0,
        "drift_time",
    )

    log.debug("attached DSP fields %s + drift_time", fields)
    return det_evt_data


def _read_dsp_field_for_channel(
    channels: ak.Array,
    rows: ak.Array,
    dsp_file: str,
    field: str,
    tab_map: dict[str, int],
) -> ak.Array:
    """Read one DSP field for all channels in a set of events using TCM indices."""
    # Validate that all rawids in the data are known
    unique_channels = np.unique(ak.flatten(channels).to_numpy())
    known_rawids = np.array(list(tab_map.values()))
    bad = unique_channels[~np.isin(unique_channels, known_rawids)]
    if len(bad) > 0:
        msg = f"channels {bad} not found in tab_map"
        raise ValueError(msg)

    # Remember the per-event structure for unflattening at the end
    counts = ak.num(rows)

    # Accumulate data and position info across channels
    data_flat = None
    tcm_rows_full = None

    for rawid in known_rawids:
        # Boolean mask: which entries in the flattened channel array match
        ch_match = channels == rawid
        idx = ak.flatten(rows[ch_match]).to_numpy()
        if len(idx) == 0:
            continue

        # lh5.read requires sorted indices
        sort_order = np.argsort(idx)
        data_ch = lh5.read(
            f"ch{rawid}/dsp/{field}", dsp_file, idx=idx[sort_order]
        ).view_as("ak")
        # Undo the sort to restore the original event order
        data_ch = data_ch[np.argsort(sort_order)]

        # Track which positions in the flat output these values belong to
        tcm_rows = np.where(ak.flatten(ch_match))[0]

        if data_flat is None:
            data_flat = data_ch
            tcm_rows_full = tcm_rows
        else:
            data_flat = ak.concatenate((data_flat, data_ch))
            tcm_rows_full = np.concatenate((tcm_rows_full, tcm_rows))

    if data_flat is None:
        msg = f"no data read for field '{field}' from {dsp_file}"
        raise ValueError(msg)

    # Re-sort into the original flattened event order and restore nesting
    data_flat = data_flat[np.argsort(tcm_rows_full)]
    return ak.unflatten(data_flat, counts)


def select_data_in_slice(
    det_evt_data: ak.Array,
    slice: Slice,
) -> ak.Array:
    """
    Filter single-detector event data to one energy-drift-time slice.

    Called after ``add_dsp_pars_to_evt``, which provides the ``drift_time``
    and energy (``geds.energy``) fields required for the slice cuts.

    The detector filter has already been applied in ``select_detector_events``,
    so this function applies only the energy and drift time bounds.

    Parameters
    ----------
    det_evt_data : ak.Array
        Single-detector event data with DSP fields attached, as returned by
        ``add_dsp_pars_to_evt``, shape ``(n_det_events,)``.
    slice : Slice
        The energy-drift-time slice to select.

    Returns
    -------
    ak.Array
        Events within the slice, shape ``(n_slice_events,)``.
    """
    mask = (
        ak.all(det_evt_data.geds.energy >= slice.energy_range[0], axis=-1)
        & ak.all(det_evt_data.geds.energy <= slice.energy_range[1], axis=-1)
        & ak.all(det_evt_data.drift_time >= slice.drift_time_range[0], axis=-1)
        & ak.all(det_evt_data.drift_time <= slice.drift_time_range[1], axis=-1)
    )
    slice_data = det_evt_data[mask]

    log.debug("%s: %d events", slice, len(slice_data))
    return slice_data


# ---------------------------------------------------------------------------
# Waveform extraction (one slice)
# ---------------------------------------------------------------------------


def get_charge_and_current_wfs_for_slice(  ### Check entire function
    raw_file: Path | str,
    lh5_group: str,
    indices: list[int],
    dsp_config: Path | str,
    charge_output: str = "wf_blsub",  ### What parameter should I use?? This? Or wf_av? Others?
    current_output: str = "curr_av",  ### What parameter should I use??
    align: str = "tp_aoe_max",
) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray | None, np.ndarray | None
]:
    """
    Extract aligned charge and current waveforms for a set of slice events, using the production DSP processing chain via ``WaveformBrowser``.

    A single ``WaveformBrowser`` instance is created for the file, with both
    ``charge_output`` and ``current_output`` requested as ``lines``. The full
    set of raw-tier row indices is passed as ``entry_list`` at construction
    time so that the browser maps them correctly to its internal buffer.
    Events are then iterated by local position (0, 1, 2, …) via
    ``find_entry`` and the processed waveforms are extracted from the
    browser's internal ``Line2D`` objects and stacked into arrays.

    The ``dsp_config`` and ``lh5_group`` required by this function can be
    obtained from ``legendsimflow.hpge_pars.lookup_currmod_fit_inputs``, which
    resolves both from the L200 data production directory given a run ID and
    detector name. The ``indices`` correspond to the raw-tier row indices for
    the slice events, obtained from ``det_evt_data.geds.hit_idx`` after
    ``select_data_in_slice``.

    Parameters
    ----------
    raw_file : Path or str
        Path to the raw-tier LH5 file containing the waveforms.
    lh5_group : str
        HDF5 group containing the waveform table, e.g. ``"ch1084803/raw"``.
    indices : list[int]
        Raw-tier row indices of the slice events to extract. These must all
        belong to ``raw_file``. Passed as ``entry_list`` to
        ``WaveformBrowser`` so that the browser can resolve them against its
        internal buffer; iteration then uses local positions 0, 1, 2, …
    dsp_config : Path or str
        Path to the production DSP configuration JSON file defining the
        processing chain. Returned as ``dsp_cfg_file`` by
        ``legendsimflow.hpge_pars.lookup_currmod_fit_inputs``.
    charge_output : str, optional
        Name of the DSP output corresponding to the baseline-subtracted charge
        waveform. Default ``"wf_blsub"``.
    current_output : str, optional
        Name of the DSP output corresponding to the MWA'd current waveform.
        Default ``"curr_av"``.
    align : str, optional
        Name of the DSP parameter used to align waveforms on the time axis.
        Default ``"tp_aoe_max"``.

    Returns
    -------
    charge_times : np.ndarray
        1D array of shape ``(n_charge_samples,)`` in ns. Time axis for the
        charge waveforms, aligned so that ``tp_aoe_max = 0``.
    current_times : np.ndarray
        1D array of shape ``(n_current_samples,)`` in ns. Time axis for the
        current waveforms, aligned so that ``tp_aoe_max = 0``.
    charge_wfs : np.ndarray
        2D array of shape ``(n_valid_events, n_charge_samples)``.
        Baseline-subtracted charge waveforms, normalised by ``cuspEmax`` so
        the plateau equals 1, aligned to ``tp_aoe_max``.
    current_wfs : np.ndarray
        2D array of shape ``(n_valid_events, n_current_samples)``. MWA'd
        current waveforms aligned to ``tp_aoe_max``.
    bl_std : np.ndarray or None
        1D array of shape ``(n_valid_events,)`` with the baseline standard
        deviation (in ADC units) for each event, as computed by the DSP chain.
        ``None`` if the DSP chain does not provide ``bl_std``.
        Used downstream by ``compute_chi2_vs_superpulse`` to estimate the
        per-event noise: ``sigma = bl_std / cuspEmax``.
    cuspEmax : np.ndarray or None
        1D array of shape ``(n_valid_events,)`` with the energy estimator
        (in ADC units) for each event. ``None`` if the DSP chain does not
        provide ``cuspEmax``. Needed alongside ``bl_std`` to compute the
        normalised noise estimate.

    Notes
    -----
    Both ``bl_std`` and ``cuspEmax`` are retrieved via the browser's
    ``legend`` mechanism as scalar DSP outputs. Since ``WaveformBrowser``
    wraps them as ``pint.Quantity`` objects, their numeric values are
    extracted via ``.magnitude``.

    Events for which the browser returns NaN-valued waveforms are dropped.
    The returned arrays contain only valid events.

    Raises
    ------
    ValueError
        If ``indices`` is empty or no valid waveforms are found across all
        ``indices``.
    """
    # HACK: deferred import to avoid pint unit registry conflict at module level
    from dspeed.vis import WaveformBrowser  # noqa: PLC0415

    if not indices:
        msg = "indices list is empty"
        raise ValueError(msg)

    browser = WaveformBrowser(
        str(raw_file),
        lh5_group,
        dsp_config=str(dsp_config),
        entry_list=indices,
        lines=[charge_output, current_output],
        norm="cuspEmax",
        align=align,
        legend=["{bl_std}", "{cuspEmax}"],
    )

    charge_list: list[np.ndarray] = []
    current_list: list[np.ndarray] = []
    bl_std_list: list[float] = []
    cuspEmax_list: list[float] = []
    charge_times = None
    current_times = None

    for local_idx in range(len(indices)):
        browser.find_entry(local_idx, append=False)

        charge_lines = browser.lines[charge_output]
        current_lines = browser.lines[current_output]

        if len(charge_lines) == 0 or len(current_lines) == 0:
            log.debug(
                "entry %d (row %d): no waveform returned, skipping",
                local_idx,
                indices[local_idx],
            )
            continue

        charge_y = charge_lines[0].get_ydata()
        current_y = current_lines[0].get_ydata()

        # Extract the specific X-axis for this event
        this_charge_x = charge_lines[0].get_xdata()
        this_current_x = current_lines[0].get_xdata()

        if np.any(np.isnan(charge_y)) or np.any(np.isnan(current_y)):
            log.debug(
                "entry %d (row %d): NaN in waveform, skipping",
                local_idx,
                indices[local_idx],
            )
            continue

        ### Charge: Trust WaveformBrowser ??
        this_charge_x = charge_lines[0].get_xdata()

        ### Current: Manually reconstruct the X-axis to force alignment ??
        # We know the sampling rate is exactly 1 ns. We force the peak to sit at 0 ns.
        peak_idx = np.nanargmax(current_y)
        this_current_x = np.arange(len(current_y), dtype=float) - peak_idx

        # Capture the master grids on the first valid event
        if charge_times is None:
            charge_times = this_charge_x.copy()
        if current_times is None:
            # Create a master grid for the current (-2000 ns to +2000 ns) to ensure all slices across the entire detector share the exact same X-axis.
            current_times = np.arange(-2000.0, 2000.0, 1.0)

        # Interpolate this event's Y-values onto the master grids
        charge_y_aligned = np.interp(
            charge_times, this_charge_x, charge_y, left=np.nan, right=np.nan
        )
        current_y_aligned = np.interp(
            current_times, this_current_x, current_y, left=np.nan, right=np.nan
        )

        charge_list.append(charge_y_aligned)
        current_list.append(current_y_aligned)

        bl_std_vals = browser.legend_vals.get("bl_std", [])
        cuspEmax_vals = browser.legend_vals.get("cuspEmax", [])
        if bl_std_vals:
            bl_std_list.append(float(bl_std_vals[0].magnitude))
        if cuspEmax_vals:
            cuspEmax_list.append(float(cuspEmax_vals[0].magnitude))

    if not charge_list:
        msg = f"no valid waveforms found in {raw_file} for {len(indices)} indices"
        raise ValueError(msg)

    log.debug(
        "extracted %d/%d valid waveforms from %s",
        len(charge_list),
        len(indices),
        Path(raw_file).name,
    )

    charge_wfs = np.array(charge_list)
    current_wfs = np.array(current_list)
    bl_std = np.array(bl_std_list) if bl_std_list else None
    cuspEmax = np.array(cuspEmax_list) if cuspEmax_list else None

    return charge_times, current_times, charge_wfs, current_wfs, bl_std, cuspEmax


# ---------------------------------------------------------------------------
# Superpulse computation (one slice)
# ---------------------------------------------------------------------------


def compute_superpulse(
    charge_times: np.ndarray,
    current_times: np.ndarray,
    charge_wfs: np.ndarray,
    current_wfs: np.ndarray,
    slice: Slice,
    detector: str,
    n_events_preliminary: int,
) -> Superpulse:
    """
    Compute the average charge and current waveforms for one slice.

    Called twice per slice: once on all preprocessed waveforms (preliminary
    superpulse), and once on the golden subset after ``apply_chi2_cut``
    (final superpulse).

    Parameters
    ----------
    charge_times : np.ndarray
        Time axis for the charge waveforms in ns, shape ``(n_charge_samples,)``.
    current_times : np.ndarray
        Time axis for the current waveforms in ns, shape ``(n_current_samples,)``.
    charge_wfs : np.ndarray
        2D array of shape ``(n_events, n_charge_samples)``.
    current_wfs : np.ndarray
        2D array of shape ``(n_events, n_current_samples)``.
    slice : Slice
        The slice this superpulse belongs to.
    detector : str
        Detector name.
    n_events_preliminary : int
        Total number of events before the chi2 cut. Pass
        ``charge_wfs.shape[0]`` on the first call. On the second call, pass
        the value from the preliminary superpulse to preserve provenance.

    Returns
    -------
    Superpulse
        With ``charge_wf = np.nanmean(charge_wfs, axis=0)``,
        ``current_wf = np.nanmean(current_wfs, axis=0)``,
        ``n_events_final = charge_wfs.shape[0]``.
    """
    return Superpulse(
        charge_wf=np.nanmean(charge_wfs, axis=0),
        current_wf=np.nanmean(current_wfs, axis=0),
        charge_time_axis=charge_times,
        current_time_axis=current_times,
        slice=slice,
        detector=detector,
        n_events_preliminary=n_events_preliminary,
        n_events_final=charge_wfs.shape[0],
    )


def compute_chi2_vs_superpulse(
    charge_wfs: np.ndarray,
    superpulse: Superpulse,
    baseline_region_mask: np.ndarray | None = None,
    bl_std: np.ndarray | None = None,
    cuspEmax: np.ndarray | None = None,
) -> np.ndarray:
    """
    Compute the reduced chi-squared of each charge waveform vs the superpulse.

    Two modes for estimating the per-event noise sigma:

    - **Preferred**: if ``bl_std`` and ``cuspEmax`` are provided (from the DSP
      chain via ``get_charge_and_current_wfs_for_slice``), the normalised noise
      is ``sigma_wf = bl_std / cuspEmax`` for each event, and the superpulse
      sigma is the mean of these values.
    - **Fallback**: if either is ``None``, sigma is estimated from the standard
      deviation of each waveform in the baseline region (samples before the
      signal onset).

    Parameters
    ----------
    charge_wfs : np.ndarray
        2D array of shape ``(n_events, n_samples)``.
    superpulse : Superpulse
        Preliminary superpulse from ``compute_superpulse``.
    baseline_region_mask : np.ndarray or None, optional
        Boolean mask of shape ``(n_samples,)`` selecting the baseline region.
        If ``None``, derived from
        ``superpulse.time_axis < -superpulse.slice.drift_time_range[1]``.
        Only used in fallback mode.
    bl_std : np.ndarray or None, optional
        Per-event baseline standard deviation in ADC units, shape
        ``(n_events,)``. From ``get_charge_and_current_wfs_for_slice``.
    cuspEmax : np.ndarray or None, optional
        Per-event energy estimator in ADC units, shape ``(n_events,)``.
        From ``get_charge_and_current_wfs_for_slice``.

    Returns
    -------
    np.ndarray
        1D array of shape ``(n_events,)`` with reduced chi-squared values.

    Notes
    -----
    The chi-squared is computed as::

        chi2 = nansum((wf - superpulse.charge_wf)**2
                      / (sigma_wf**2 + sigma_sp**2))
        reduced_chi2 = chi2 / n_dof

    where ``n_dof = n_samples``.
    """
    _, n_samples = charge_wfs.shape
    sp_wf = superpulse.charge_wf

    use_dsp_noise = bl_std is not None and cuspEmax is not None

    if use_dsp_noise:
        # Preferred: normalised noise from DSP chain
        sigma_wfs = bl_std / cuspEmax  # (n_events,)
        sigma_sp = np.mean(sigma_wfs)
    else:  ### Safe enough to remove the fallback, and rely only on the value from proc chain??
        # Fallback: estimate from baseline region
        if baseline_region_mask is None:
            baseline_region_mask = (
                superpulse.charge_time_axis < -superpulse.slice.drift_time_range[1]
            )
        if baseline_region_mask.sum() < 10:
            log.warning(
                "only %d baseline samples available for noise estimation; "
                "consider providing bl_std and cuspEmax from the DSP chain",
                baseline_region_mask.sum(),
            )
        if baseline_region_mask.sum() == 0:
            msg = (
                "no baseline samples available for noise estimation; "
                "provide bl_std and cuspEmax from the DSP chain, or pass an explicit baseline_region_mask"
            )
            raise ValueError(msg)
        sigma_wfs = np.std(charge_wfs[:, baseline_region_mask], axis=1)  # (n_events,)
        sigma_sp = np.std(sp_wf[baseline_region_mask])

    # Vectorised chi2: residuals shape (n_events, n_samples), sigma_wfs broadcast from (n_events, 1)
    residuals_sq = (charge_wfs - sp_wf[np.newaxis, :]) ** 2
    variance = sigma_wfs[:, np.newaxis] ** 2 + sigma_sp**2
    chi2 = np.nansum(residuals_sq / variance, axis=1)

    return chi2 / n_samples


def apply_chi2_cut(
    charge_wfs: np.ndarray,
    current_wfs: np.ndarray,
    chi2_values: np.ndarray,
    threshold: float = 3.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Retain only waveforms with reduced chi-squared below ``threshold``.

    Parameters
    ----------
    charge_wfs : np.ndarray
        2D array of shape ``(n_events, n_samples)``.
    current_wfs : np.ndarray
        2D array of shape ``(n_events, n_samples)``.
    chi2_values : np.ndarray
        1D array of reduced chi-squared values from
        ``compute_chi2_vs_superpulse``.
    threshold : float, optional
        Reduced chi-squared cut value. Default 3.0.

    Returns
    -------
    golden_charge_wfs : np.ndarray
        Shape ``(n_golden, n_samples)``.
    golden_current_wfs : np.ndarray
        Shape ``(n_golden, n_samples)``.
    golden_indices : np.ndarray
        1D integer array of shape ``(n_golden,)`` with the indices of
        surviving events in the input arrays. Retained for traceability
        back to the original ``slice_data_with_wfs``.
    """
    golden_indices = np.where(chi2_values < threshold)[0]
    log.debug(
        "chi2 cut (< %.1f): %d / %d waveforms pass",
        threshold,
        len(golden_indices),
        len(chi2_values),
    )
    return charge_wfs[golden_indices], current_wfs[golden_indices], golden_indices


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def build_superpulse_for_slice(
    slice: Slice,
    detector: str,
    tab_map: dict[str, int],
    evt_files: list[str],
    dsp_files: list[str],
    raw_files: list[str],
    dsp_config: Path | str,
    dsp_fields: list[str],
    n_target_wfs: int = 100,
    chi2_threshold: float = 3.0,
    charge_output: str = "wf_av",
    current_output: str = "curr_av",
    align: str = "tp_aoe_max",
) -> Superpulse:
    """
    Build the final superpulse for one slice using an adaptive two-phase loop.

    This is the main entry point for superpulse construction. It manages the
    file loop and the two-phase adaptive strategy:

    **Phase 1 - accumulate:** read files one at a time, running the full
    dataset preparation pipeline (steps 1-6) on each. Accumulate the
    resulting charge and current waveform arrays. Continue until the
    accumulated event count reaches ``n_target_wfs``.

    **Phase 2 - build and check:** run the full superpulse construction
    pipeline (steps 7-10) on the accumulated data. If the number of golden
    waveforms after the chi2 cut is below ``n_target_wfs``, return to
    Phase 1 and continue reading from the next unread file, adding to the
    existing accumulation. Repeat until either:

    - ``n_golden >= n_target_wfs``: success, return the final superpulse.
    - All files are exhausted: log a warning and return the best superpulse
      available, or raise if no waveforms were found at all.

    The ``lh5_group`` (HDF5 group containing the waveform table, e.g. ``"ch1084803/raw"``)
    and ``dsp_config`` parameters are most conveniently
    obtained from ``legendsimflow.hpge_pars.lookup_currmod_fit_inputs``,
    which resolves both from the L200 data production directory given a
    run ID and detector name.

    Parameters
    ----------
    slice : Slice
        The energy-drift-time slice to process.
    detector : str
        Name of the target detector, e.g. ``"V03422A"``.
    tab_map : dict[str, int]
        Mapping from detector name to rawid.
    evt_files : list[str]
        Sorted list of all available EVT-tier LH5 file paths.
    dsp_files : list[str]
        Sorted list of all available DSP-tier LH5 file paths, same order as
        ``evt_files``.
    raw_files : list[str]
        Sorted list of all available raw-tier LH5 file paths, same order as
        ``evt_files``.
    dsp_config : Path or str
        Path to the production DSP configuration JSON file. Returned as
        ``dsp_cfg_file`` by
        ``legendsimflow.hpge_pars.lookup_currmod_fit_inputs``.
    dsp_fields : list[str]
        DSP fields to attach to events via ``add_dsp_pars_to_evt``.
        ``tp_aoe_max`` must be included.
    n_target_wfs : int, optional
        Target number of golden waveforms (after chi2 cut) for the final
        superpulse. Default 100.
    chi2_threshold : float, optional
        Reduced chi-squared threshold for the self-similarity cut. Default 3.0.
    charge_output : str, optional
        DSP output name for the charge waveform. Default ``"wf_blsub"``.
        Passed to ``get_charge_and_current_wfs_for_slice``.
    current_output : str, optional
        DSP output name for the current waveform. Default ``"curr_av"``.
        Passed to ``get_charge_and_current_wfs_for_slice``.
    align : str, optional
        DSP parameter used for waveform alignment. Default ``"tp_aoe_max"``.
        Passed to ``get_charge_and_current_wfs_for_slice``.

    Returns
    -------
    Superpulse
        Final superpulse with ``n_events_final >= n_target_wfs`` if enough
        data was available, otherwise the best achievable superpulse from all
        available files.

    Raises
    ------
    ValueError
        If no waveforms are found in the slice across all files.
    RuntimeError
        If ``evt_files``, ``dsp_files``, and ``raw_files`` have different
        lengths.

    Notes
    -----
    File triplets ``(evt_files[i], dsp_files[i], raw_files[i])`` must
    correspond to the same run segment. The caller is responsible for
    ensuring this ordering.

    Memory usage is proportional to the number of accumulated events times
    the waveform length. Since waveforms are extracted already upsampled and
    processed by the DSP chain, no further upsampling is performed in this
    function.
    """
    # Input validation
    if len(evt_files) != len(dsp_files) or len(evt_files) != len(raw_files):
        msg = (
            f"file lists must have equal length: "
            f"evt={len(evt_files)}, dsp={len(dsp_files)}, raw={len(raw_files)}"
        )
        raise RuntimeError(msg)

    # Construct the raw-tier LH5 group from the rawid
    rawid = tab_map[detector]
    lh5_group = f"ch{rawid}/raw"

    # Accumulators
    all_charge_wfs: list[np.ndarray] = []
    all_current_wfs: list[np.ndarray] = []
    all_bl_std: list[np.ndarray] = []
    all_cuspEmax: list[np.ndarray] = []
    charge_common_times: np.ndarray | None = None
    current_common_times: np.ndarray | None = None

    # Overshoot factor: accumulate more raw waveforms than needed,
    # since the chi2 cut will reject some fraction
    _ACCUMULATION_FACTOR = 1.5  ### I can increase it to avoid looping back to phase 1, but it will increase execution time and memory usage

    file_idx = 0  # tracks which file to read next

    while file_idx < len(evt_files):
        # ---------------------------------------------------------------
        # Phase 1: accumulate waveforms from unread files
        # ---------------------------------------------------------------
        while file_idx < len(evt_files):
            log.info(
                "%s | %s | file %d/%d",
                detector,
                slice,
                file_idx + 1,
                len(evt_files),
            )

            # read EVT data for one file
            evt_data = read_evt_data([evt_files[file_idx]])

            # quality + PSD cuts
            evt_data = perform_data_selection(evt_data)

            # filter to target detector
            det_data = select_detector_events(evt_data, detector, tab_map)

            if len(det_data) == 0:
                log.debug("file %d: no events for %s after cuts", file_idx, detector)
                file_idx += 1
                continue

            # attach DSP fields + drift_time
            det_data = add_dsp_pars_to_evt(
                det_data,
                dsp_files[file_idx],
                detector,
                tab_map,
                dsp_fields,
            )

            # select events in the slice
            slice_data = select_data_in_slice(det_data, slice)

            if len(slice_data) == 0:
                log.debug("file %d: no events in %s", file_idx, slice)
                file_idx += 1
                continue

            # extract waveforms via WaveformBrowser
            indices = ak.flatten(slice_data.geds.hit_idx).to_list()

            try:
                (
                    charge_times,
                    current_times,
                    charge_wfs,
                    current_wfs,
                    bl_std,
                    cuspEmax,
                ) = get_charge_and_current_wfs_for_slice(
                    raw_file=raw_files[file_idx],
                    lh5_group=lh5_group,
                    indices=indices,
                    dsp_config=dsp_config,
                    charge_output=charge_output,
                    current_output=current_output,
                    align=align,
                )
            except ValueError:
                log.debug("file %d: no valid waveforms extracted", file_idx)
                file_idx += 1
                continue

            # Store the common time axes (same for all files after alignment)
            if charge_common_times is None:
                charge_common_times = charge_times
            if current_common_times is None:
                current_common_times = current_times

            all_charge_wfs.append(charge_wfs)
            all_current_wfs.append(current_wfs)
            if bl_std is not None:
                all_bl_std.append(bl_std)
            if cuspEmax is not None:
                all_cuspEmax.append(cuspEmax)

            file_idx += 1

            # Check if we have enough raw waveforms to attempt the chi2 cut
            n_accumulated = sum(c.shape[0] for c in all_charge_wfs)
            if n_accumulated >= _ACCUMULATION_FACTOR * n_target_wfs:
                break

        # ---------------------------------------------------------------
        # Phase 2: build superpulse and apply chi2 cut
        # ---------------------------------------------------------------
        n_accumulated = sum(c.shape[0] for c in all_charge_wfs)
        if n_accumulated == 0:
            # No waveforms found yet; if files remain, keep going
            continue

        # Stack all accumulated waveforms
        charge_stack = np.concatenate(all_charge_wfs, axis=0)
        current_stack = np.concatenate(all_current_wfs, axis=0)
        bl_std_stack = np.concatenate(all_bl_std, axis=0) if all_bl_std else None
        cuspEmax_stack = np.concatenate(all_cuspEmax, axis=0) if all_cuspEmax else None

        n_preliminary = charge_stack.shape[0]

        # preliminary superpulse
        prelim_sp = compute_superpulse(
            charge_common_times,
            current_common_times,
            charge_stack,
            current_stack,
            slice,
            detector,
            n_preliminary,
        )

        # chi2 of each waveform vs preliminary superpulse
        chi2_values = compute_chi2_vs_superpulse(
            charge_stack,
            prelim_sp,
            bl_std=bl_std_stack,
            cuspEmax=cuspEmax_stack,
        )

        # apply chi2 cut
        golden_charge, golden_current, _ = apply_chi2_cut(
            charge_stack,
            current_stack,
            chi2_values,
            threshold=chi2_threshold,
        )
        n_golden = golden_charge.shape[0]

        log.info(
            "%s | %s | golden: %d / %d (target: %d)",
            detector,
            slice,
            n_golden,
            n_preliminary,
            n_target_wfs,
        )

        if n_golden >= n_target_wfs or file_idx >= len(evt_files):
            # Either we have enough, or we've exhausted all files.
            # In both cases: build the final superpulse from golden wfs.
            if n_golden == 0:
                log.warning(
                    "%s | %s | no waveforms survived chi2 cut; "
                    "returning preliminary superpulse",
                    detector,
                    slice,
                )
                return prelim_sp

            if n_golden < n_target_wfs:
                log.warning(
                    "%s | %s | all files exhausted: %d golden waveforms "
                    "(target was %d)",
                    detector,
                    slice,
                    n_golden,
                    n_target_wfs,
                )

            return compute_superpulse(
                charge_common_times,
                current_common_times,
                golden_charge,
                golden_current,
                slice,
                detector,
                n_preliminary,
            )

        # Not enough golden waveforms and files remain - loop back to Phase 1
        log.info(
            "%s | %s | need more waveforms, continuing to next file...",
            detector,
            slice,
        )

    # Only reachable if every file yielded zero waveforms for this slice
    msg = f"no waveforms found for {detector} in {slice} across all files"
    raise ValueError(msg)


# ---------------------------------------------------------------------------
# Output: all slices, one detector
# ---------------------------------------------------------------------------


def write_superpulses_to_lh5(
    superpulses: dict[Slice, Superpulse],
    output_path: str,
    detector: str,
) -> None:
    """
    Write all per-slice superpulses for one detector to a single LH5 file.

    Iterates over ``superpulses``, calls ``Superpulse.to_lgdo()`` for each,
    and writes the result as a named group inside the file.

    The LH5 file structure is::

        {detector}/
            dt_{lo}_{hi}_ns/
                charge_wf             [Array, n_charge_samples]
                current_wf            [Array, n_current_samples]
                charge_time_axis      [Array, n_charge_samples, units=ns]
                current_time_axis     [Array, n_current_samples, units=ns]
                dt_center             [Scalar, ns]
                dt_lo                 [Scalar, ns]
                dt_hi                 [Scalar, ns]
                e_lo                  [Scalar, keV]
                e_hi                  [Scalar, keV]
                detector              [Scalar, str]
                n_events_preliminary  [Scalar]
                n_events_final        [Scalar]

    Parameters
    ----------
    superpulses : dict[Slice, Superpulse]
        Dictionary mapping each slice to its final superpulse, as accumulated
        in the per-slice processing loop.
    output_path : str
        Path to the output LH5 file, e.g.
        ``"output/V03422A_superpulses.lh5"``.
    detector : str
        Detector name, used as the top-level group name in the LH5 file.

    Notes
    -----
    Uses ``lgdo.lh5.write`` with ``wo_mode="write_safe"`` to avoid silently
    overwriting existing files.
    """
    for sl, sp in superpulses.items():
        # Build group name: {detector}/dt_{lo}_{hi}_ns
        dt_lo = int(sl.drift_time_range[0])
        dt_hi = int(sl.drift_time_range[1])
        group = f"{detector}/dt_{dt_lo}_{dt_hi}_ns"

        lgdo_struct = sp.to_lgdo()
        lh5.write(lgdo_struct, group, output_path, wo_mode="write_safe")

        log.debug("wrote %s to %s", group, output_path)

    log.info(
        "wrote %d slices for %s to %s",
        len(superpulses),
        detector,
        output_path,
    )
