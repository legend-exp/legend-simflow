"""
This module implements the dataset preparation and superpulse construction
for SSC data.

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
Steps 1-7 are called once per file inside ``build_superpulse_for_pixel``,
which accumulates results across files until enough waveforms are available.

1. ``read_evt_data``               - read EVT tier for one file (no DSP fields)
2. ``perform_data_selection``      - quality + PSD cuts on EVT fields only
3. ``select_detector_events``      - filter to one detector
4. ``add_dsp_pars_to_evt``         - attach DSP fields for one detector only
5. ``select_data_in_pixel``        - filter to one energy-drift-time pixel
6. ``read_wfs_for_pixel``          - read raw waveforms for pixel events only
7. ``merge_wfs_into_evt``          - attach waveforms to event array

Steps 8-12 run once, after enough waveforms have been accumulated:

8.  ``preprocess_waveforms``        - normalize, upsample, differentiate, MWA,
                                      recalculate tp_aoe_max, align, trim
9.  ``compute_superpulse``          - preliminary average (charge + current)
10. ``compute_chi2_vs_superpulse``  - reduced chi2 of each wf vs superpulse
11. ``apply_chi2_cut``              - retain golden waveforms only
12. ``compute_superpulse``          - final average on golden waveforms

If after step 11 the number of golden waveforms is below ``n_target_wfs``,
the orchestrator returns to step 1 and continues reading files from where
it left off, re-running steps 8-12 on the enlarged accumulated set, until
either the target is reached or all files are exhausted.

Step 13 runs once per detector, after all pixels are done:

13. ``write_superpulses_to_lh5``   - write all pixels for one detector to LH5

Dependencies
------------
- lgdo  
- awkward
- numpy 
- legendmeta
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import awkward as ak
import numpy as np
from lgdo import lh5

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Pixel:
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
            f"Pixel(E=[{self.energy_range[0]}, {self.energy_range[1]}] keV, "
            f"dt=[{self.drift_time_range[0]}, {self.drift_time_range[1]}] ns)"
        )


class Superpulse:
    """
    Holds the average charge and current waveforms for one pixel.

    One ``Superpulse`` instance is produced per pixel per detector, after the
    full preprocessing and chi2 self-similarity cut. It carries both waveforms
    together with the metadata needed to identify the pixel, write to LH5, and
    perform the subsequent electronics parameter optimisation.

    Parameters
    ----------
    charge_wf : np.ndarray
        Average normalised charge waveform, shape ``(n_samples,)``.
        Amplitude is dimensionless (ADC / cuspEmax, normalised to 1 at
        plateau).
    current_wf : np.ndarray
        Average current waveform after MWA, shape ``(n_samples,)``.
        Units: ``(ADC / cuspEmax) / ns * dt_data``, matching the convention
        in ``psl.py``.
    time_axis : np.ndarray
        Common time axis in ns, shape ``(n_samples,)``, aligned so that
        ``tp_aoe_max = 0``.
    pixel : Pixel
        The energy-drift-time pixel this superpulse represents.
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
            time_axis=t,
            pixel=Pixel((1500., 2000.), (900., 1100.)),
            detector="V03422A",
            n_events_preliminary=120,
            n_events_final=98,
        )
        sp.charge_wf            # np.ndarray
        sp.current_wf           # np.ndarray
        sp.pixel.drift_time_center  # 1000.0 ns

    Use as dict value (Pixel is hashable)::

        superpulses: dict[Pixel, Superpulse] = {}
        superpulses[sp.pixel] = sp
    """

    def __init__(
        self,
        charge_wf: np.ndarray,
        current_wf: np.ndarray,
        time_axis: np.ndarray,
        pixel: Pixel,
        detector: str,
        n_events_preliminary: int,
        n_events_final: int,
    ) -> None:
        if len(charge_wf) != len(time_axis) or len(current_wf) != len(time_axis):
            msg = (
                f"charge_wf ({len(charge_wf)}), current_wf ({len(current_wf)}), "
                f"and time_axis ({len(time_axis)}) must have the same length."
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
        self.time_axis = time_axis
        self.pixel = pixel
        self.detector = detector
        self.n_events_preliminary = n_events_preliminary
        self.n_events_final = n_events_final

    def __repr__(self) -> str:
        return (
            f"Superpulse("
            f"detector={self.detector!r}, "
            f"pixel={self.pixel}, "
            f"n_events={self.n_events_final}/{self.n_events_preliminary}, "
            f"n_samples={len(self.time_axis)})"
        )

    def to_lgdo(self):
        """
        Serialise to an ``lgdo.Struct`` ready for writing to LH5.

        Returns
        -------
        lgdo.Struct
            Struct with the following fields:

            - ``charge_wf``            : ``lgdo.Array``, shape ``(n_samples,)``
            - ``current_wf``           : ``lgdo.Array``, shape ``(n_samples,)``
            - ``time_axis``            : ``lgdo.Array``, shape ``(n_samples,)``,
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
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Dataset preparation (individual steps, one file at a time)
# ---------------------------------------------------------------------------


def read_evt_data(
    evt_files: list[str],
) -> ak.Array:
    """
    Read EVT-tier LH5 files into a single concatenated awkward array.

    No DSP fields are attached here. This is intentionally a lightweight
    read so that ``perform_data_selection`` and ``select_detector_events``
    can reduce the array size before the DSP read in ``add_dsp_pars_to_evt``.

    In normal use this is called with a single-element list by the orchestrator
    ``build_superpulse_for_pixel``, which manages the file loop externally.
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
    raise NotImplementedError


def perform_data_selection(
    evt_data: ak.Array,
    aoe_low_threshold: float = -3.0,
    aoe_high_threshold: float = 3.0,
) -> ak.Array:
    """
    Apply quality and PSD cuts to retain only valid single-site events.

    Two sequential sets of cuts are applied:

    1. **Preliminary cuts**: forced triggers, pulser events, muon coincidences
       (online and offline), bad channel quality, failing QC, multiplicity > 1,
       and insufficient LAr light (``spms.energy_sum <= 10``).

    2. **PSD cuts**: A/E lower and upper bounds on ``geds.psd.low_aoe.value``
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
    raise NotImplementedError


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
    raise NotImplementedError


def add_dsp_pars_to_evt(
    det_evt_data: ak.Array,
    dsp_files: list[str],
    detector: str,
    tab_map: dict[str, int],
    fields: list[str],
) -> ak.Array:
    """
    Attach per-channel DSP parameters to a single-detector event array.

    Called after ``select_detector_events``, so ``det_evt_data`` already
    contains only events for one detector. The DSP read therefore involves
    only one channel: one ``lh5.read`` call per field per file.

    ``drift_time = tp_aoe_max - spms.first_t0`` is also computed and added
    as a top-level field, since it requires both a DSP field (``tp_aoe_max``)
    and an EVT field (``spms.first_t0``).

    In normal use ``dsp_files`` is a single-element list corresponding to the
    same run segment as the EVT file passed to ``read_evt_data``.

    Parameters
    ----------
    det_evt_data : ak.Array
        Single-detector event data as returned by ``select_detector_events``,
        shape ``(n_det_events,)``.
    dsp_files : list[str]
        DSP-tier LH5 file path(s), in the same order as the EVT files used
        to produce ``det_evt_data``.
    detector : str
        Name of the target detector, e.g. ``"V03422A"``.
    tab_map : dict[str, int]
        Mapping from detector name to rawid. Only the entry for ``detector``
        is used.
    fields : list[str]
        DSP fields to attach, e.g. ``["tp_aoe_max", "cuspEmax", "tp_0_est"]``.
        ``tp_aoe_max`` must be included for ``drift_time`` to be computed.

    Returns
    -------
    ak.Array
        ``det_evt_data`` with the requested DSP fields and ``drift_time``
        added as top-level fields. Each field has shape ``(1,)`` per event,
        consistent with the multiplicity == 1 guarantee.

    Notes
    -----
    The TCM mapping uses ``det_evt_data.geds.rawid`` and
    ``det_evt_data.geds.hit_idx`` to identify which DSP rows correspond to
    each event, via ``_read_dsp_field_for_channel``.

    Raises
    ------
    ValueError
        If ``tp_aoe_max`` is not in ``fields`` (required for drift_time).
    """
    raise NotImplementedError


def _read_dsp_field_for_channel(
    channels: ak.Array,
    rows: ak.Array,
    dsp_file: str,
    field: str,
    tab_map: dict[str, int],
) -> ak.Array:
    """
    Read one DSP field for all channels in a set of events, using TCM indices.

    Low-level helper for ``add_dsp_pars_to_evt``. For each channel in
    ``tab_map``, identifies the corresponding TCM rows, reads only those rows
    from the DSP file via ``lgdo.lh5.read(..., idx=sorted_idx)``, then
    re-orders results to match the original event ordering.

    When called from a single-detector context (as intended in this pipeline),
    ``tab_map`` has one entry and the channel loop executes exactly once.

    Parameters
    ----------
    channels : ak.Array
        Shape ``(n_events, n_hits_per_event)``. ``evt_data.geds.rawid``.
    rows : ak.Array
        Shape ``(n_events, n_hits_per_event)``. ``evt_data.geds.hit_idx``.
    dsp_file : str
        Path to the DSP-tier LH5 file.
    field : str
        DSP field name, e.g. ``"tp_aoe_max"``.
    tab_map : dict[str, int]
        Mapping from detector name to rawid.

    Returns
    -------
    ak.Array
        Shape ``(n_events, n_hits_per_event)``, values of ``field`` in
        original event order.

    Notes
    -----
    ``lh5.read()`` requires a sorted ``idx`` array. The default
    ``use_h5idx=False`` reads the full column and then indexes in memory:
    faster than random HDF5 row access for most sizes, but uses more memory.
    For very sparse reads this trade-off should be re-evaluated.

    Raises
    ------
    ValueError
        If any rawid in ``channels`` is absent from ``tab_map``.
    """
    raise NotImplementedError


def select_data_in_pixel(
    det_evt_data: ak.Array,
    pixel: Pixel,
) -> ak.Array:
    """
    Filter single-detector event data to one energy-drift-time pixel.

    Called after ``add_dsp_pars_to_evt``, which provides the ``drift_time``
    and energy (``geds.energy``) fields required for the pixel cuts.

    The detector filter has already been applied in ``select_detector_events``,
    so this function applies only the energy and drift time bounds.

    Parameters
    ----------
    det_evt_data : ak.Array
        Single-detector event data with DSP fields attached, as returned by
        ``add_dsp_pars_to_evt``, shape ``(n_det_events,)``.
    pixel : Pixel
        The energy-drift-time pixel to select.

    Returns
    -------
    ak.Array
        Events within the pixel, shape ``(n_pixel_events,)``.
    """
    raise NotImplementedError


def read_wfs_for_pixel(
    raw_files: list[str],
    pixel_evt_data: ak.Array,
    rawid: int,
) -> tuple[ak.Array, np.ndarray]:
    """
    Read raw waveforms for only the events in a given pixel.

    Uses ``trigger.timestamp`` values from ``pixel_evt_data`` to locate
    matching entries in the raw tier, then reads only those rows via
    ``lgdo.lh5.read(..., idx=sorted_idx, field_mask=[...])``.

    In normal use ``raw_files`` is a single-element list corresponding to the
    same run segment as the EVT file being processed in the current iteration
    of the orchestrator loop.

    Parameters
    ----------
    raw_files : list[str]
        Raw-tier LH5 file path(s) to search.
    pixel_evt_data : ak.Array
        Pixel event data as returned by ``select_data_in_pixel``,
        shape ``(n_pixel_events,)``.
    rawid : int
        Raw channel ID, e.g. ``1084803``. Used to build ``ch{rawid}/raw/``.

    Returns
    -------
    waveforms : ak.Array
        Shape ``(n_matched_events,)`` with fields:

        - ``timestamp`` : int64 [ns]
        - ``waveform_windowed`` : struct with ``t0`` (float [ns]),
          ``dt`` (float [ns/sample]), ``values`` (float32[n_samples])

        ``n_matched_events <= n_pixel_events`` if raw files are missing.
    valid_mask : np.ndarray
        Boolean array of shape ``(n_pixel_events,)``. Entry ``i`` is ``True``
        if ``pixel_evt_data[i]`` has a matching waveform. Use as:
        ``pixel_evt_data[valid_mask]`` aligns row-for-row with ``waveforms``.

    Notes
    -----
    Timestamp matching is the only reliable link between EVT-tier events and
    raw-tier waveforms, as the TCM does not index into the raw tier directly.

    Per file: timestamps are read first (cheap scalar read), matching row
    indices are computed, then waveform data for those rows only is read.

    The digitizer ``baseline`` field is not returned: it is imprecise and
    recomputed from waveform samples in preprocessing.

    Missing raw files are logged as warnings and skipped gracefully.

    Raises
    ------
    ValueError
        If no waveforms are found across all files.
    """
    raise NotImplementedError


def merge_wfs_into_evt(
    pixel_evt_data: ak.Array,
    waveforms: ak.Array,
    valid_mask: np.ndarray,
) -> ak.Array:
    """
    Attach waveforms to the pixel event array, retaining only matched events.

    Parameters
    ----------
    pixel_evt_data : ak.Array
        Pixel event data, shape ``(n_pixel_events,)``.
    waveforms : ak.Array
        Waveform data from ``read_wfs_for_pixel``, shape ``(n_matched,)``.
    valid_mask : np.ndarray
        Boolean mask of shape ``(n_pixel_events,)`` from ``read_wfs_for_pixel``.

    Returns
    -------
    ak.Array
        Shape ``(n_matched,)``. All fields from ``pixel_evt_data[valid_mask]``
        plus ``waveform_windowed`` (struct with ``t0``, ``dt``, ``values``).

        The digitizer ``baseline`` is not included; baseline subtraction is
        recomputed from waveform samples in preprocessing.

    Notes
    -----
    Both sides are sorted by ``trigger.timestamp`` before merging to guarantee
    row-for-row correspondence. With O(100) events per pixel this is negligible.

    Raises
    ------
    ValueError
        If ``np.sum(valid_mask) != len(waveforms)``.
    """
    raise NotImplementedError


# ---------------------------------------------------------------------------
# Waveform preprocessing (one pixel)
# ---------------------------------------------------------------------------


def preprocess_waveforms(
    pixel_data_with_wfs: ak.Array,
    bl_window: tuple[float, float],
    upsample_factor: int = 16,
    upsample_mode: str = "l",
    mwa_length: int = 48,
    mwa_num: int = 3,
    mwa_type: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Full preprocessing pipeline for one pixel: normalize, upsample,
    differentiate, MWA, recalculate tp_aoe_max, align, trim.

    The order of steps is physically fixed and must not be changed:

    1. Baseline subtraction and normalisation by ``cuspEmax``
    2. Upsample to 1 ns (factor 16 for 16 ns data)
    3. Compute derivative (current)
    4. Apply MWA to current
    5. Recalculate ``tp_aoe_max`` from MWA'd current at 1 ns resolution
    6. Align charge and current waveforms to recalculated ``tp_aoe_max``
    7. Trim all waveforms to the common overlapping time window

    Parameters
    ----------
    pixel_data_with_wfs : ak.Array
        Per-pixel event data with waveforms attached, as returned by
        ``merge_wfs_into_evt``. Each event must have:

        - ``waveform_windowed``: struct with ``t0`` (float), ``dt`` (float),
          ``values`` (float32[])
        - ``cuspEmax``: float, uncalibrated energy in ADC. Used for
          normalisation so that the charge waveform plateau equals 1.
        - ``tp_aoe_max``: float, DSP-computed current maximum time in ns.
          Used only for validation; the actual alignment uses the value
          recalculated at 1 ns resolution in step 5.

    bl_window : tuple[float, float]
        Time window ``(t_start, t_end)`` in ns in the original, non-aligned
        time axis, used to compute the per-event baseline as the mean of
        samples within the window.
    upsample_factor : int, optional
        Upsampling factor. Default 16 (16 ns raw sampling to 1 ns).
    upsample_mode : str, optional
        Interpolation mode for ``dspeed.processors.interpolating_upsampler``.
        ``'l'`` = linear, ``'h'`` = Hermite, ``'s'`` = spline. Default ``'l'``.
    mwa_length : int, optional
        Moving window length in samples at 1 ns resolution. Default 48.
        Must match the value used in ``psl.py:make_realistic_pulse_shape_lib``.
    mwa_num : int, optional
        Number of moving windows. Default 3.
        Must match the value used in ``psl.py:make_realistic_pulse_shape_lib``.
    mwa_type : int, optional
        Moving window type for ``dspeed.processors.moving_window_multi``.
        Default 0 (alternating left/right). Must match ``psl.py``.

    Returns
    -------
    common_times : np.ndarray
        1D array of shape ``(n_samples,)`` in ns. Zero corresponds to
        ``tp_aoe_max``.
    charge_wfs : np.ndarray
        2D array of shape ``(n_valid_events, n_samples)``. Normalised charge
        waveforms aligned to ``tp_aoe_max``.
    current_wfs : np.ndarray
        2D array of shape ``(n_valid_events, n_samples)``. MWA'd current
        waveforms aligned to ``tp_aoe_max``.
    valid_mask : np.ndarray
        Boolean array of shape ``(len(pixel_data_with_wfs),)``. ``True`` for
        events that survived all preprocessing sub-steps (valid ``cuspEmax``,
        successful upsampling, finite waveform values).

    Notes
    -----
    The ``geds.energy`` field (calibrated, in keV) is used for pixel selection
    upstream. Normalisation here uses ``cuspEmax`` (uncalibrated, in ADC) so
    that the charge waveform amplitude equals 1 at the plateau, independent of
    gain calibration.
    """
    raise NotImplementedError


def _subtract_baseline_and_normalise(
    wf_values: np.ndarray,
    time_axis: np.ndarray,
    cuspEmax: float,
    bl_window: tuple[float, float],
) -> tuple[np.ndarray, bool]:
    """
    Subtract baseline and normalise a single waveform by ``cuspEmax``.

    Parameters
    ----------
    wf_values : np.ndarray
        Raw waveform samples.
    time_axis : np.ndarray
        Time axis in ns corresponding to ``wf_values``, using the original
        (non-aligned) time values from ``waveform_windowed.t0`` and ``dt``.
    cuspEmax : float
        Uncalibrated energy in ADC. The waveform is divided by this value
        so that the plateau amplitude equals 1.
    bl_window : tuple[float, float]
        ``(t_start, t_end)`` in ns. Baseline is computed as the mean of
        samples with time in ``[t_start, t_end]``.

    Returns
    -------
    wf_normalised : np.ndarray
        Baseline-subtracted and energy-normalised waveform.
    valid : bool
        ``False`` if ``cuspEmax <= 0`` (invalid event, must be excluded).
    """
    raise NotImplementedError


def _upsample_waveform(
    waveform: np.ndarray,
    upsample_factor: int = 16,
    mode: str = "l",
) -> np.ndarray:
    """
    Upsample a single waveform using ``dspeed.processors.interpolating_upsampler``.

    Parameters
    ----------
    waveform : np.ndarray
        Input waveform, dtype float32.
    upsample_factor : int, optional
        Upsampling factor. Default 16.
    mode : str, optional
        Interpolation mode: ``'l'`` (linear), ``'h'`` (Hermite), ``'s'``
        (spline). Default ``'l'``.

    Returns
    -------
    np.ndarray
        Upsampled waveform of length ``len(waveform) * upsample_factor``.
    """
    raise NotImplementedError


def _compute_current(
    waveform: np.ndarray,
    dt_upsampled: float,
    dt_data: float,
) -> np.ndarray:
    """
    Compute the current waveform as a scaled discrete derivative.

    Uses ``np.diff`` with ``prepend=0``, scaled by ``dt_data / dt_upsampled``
    to preserve physical units. This mirrors the inline expression in
    ``psl.py:make_realistic_pulse_shape_lib``::

        np.diff(wfs, axis=-1, prepend=0) * (dt_data / dt)

    Parameters
    ----------
    waveform : np.ndarray
        Normalised, upsampled charge waveform.
    dt_upsampled : float
        Time step of the upsampled waveform in ns. Typically 1.0 ns.
    dt_data : float
        Original data sampling period in ns. Typically 16.0 ns.

    Returns
    -------
    np.ndarray
        Current waveform, same length as input.

    Notes
    -----
    The scaling convention here must match that used in ``psl.py``.
    If either ``dt`` value is changed, both modules must be updated
    consistently.
    """
    raise NotImplementedError


def _apply_mwa(
    current: np.ndarray,
    mwa_length: int = 48,
    mwa_num: int = 3,
    mwa_type: int = 0,
) -> np.ndarray:
    """
    Apply a multi-pass moving window average to a single current waveform,
    using ``dspeed.processors.moving_window_multi``.

    Parameters
    ----------
    current : np.ndarray
        Current waveform at 1 ns sampling, dtype float32.
    mwa_length : int, optional
        Moving window length in samples. Default 48.
    mwa_num : int, optional
        Number of passes. Default 3.
    mwa_type : int, optional
        Window type passed to ``dspeed.processors.moving_window_multi``.
        Default 0 (alternating left/right).

    Returns
    -------
    np.ndarray
        MWA-smoothed current waveform, same length as input.

    Warnings
    --------
    ``mwa_length``, ``mwa_num``, and ``mwa_type`` must be identical to the
    values in ``mw_pars`` passed to ``psl.py:make_realistic_pulse_shape_lib``.
    A mismatch will cause a systematic bias in the electronics parameter
    optimisation.

    Notes
    -----
    The ``dspeed`` import is deferred to the function body to avoid a known
    pint unit registry conflict at module level (see ``psl.py`` for the same
    pattern).
    """
    raise NotImplementedError


def _recalculate_tp_aoe_max(
    mwa_current: np.ndarray,
    time_axis: np.ndarray,
    dsp_tp_aoe_max: float,
) -> float:
    """
    Find the time of the current maximum at 1 ns resolution.

    Replaces the DSP-stored ``tp_aoe_max`` (computed at 16 ns resolution)
    with a value computed on the upsampled, MWA'd current.

    Parameters
    ----------
    mwa_current : np.ndarray
        MWA'd current waveform at 1 ns sampling.
    time_axis : np.ndarray
        Corresponding time axis in ns (not yet aligned).
    dsp_tp_aoe_max : float
        DSP-stored ``tp_aoe_max`` in ns, used only for validation logging.

    Returns
    -------
    float
        Recalculated ``tp_aoe_max`` in ns.

    Notes
    -----
    The difference ``|recalculated - dsp_tp_aoe_max|`` is logged at DEBUG
    level for validation. A large systematic difference (>> 16 ns) would
    indicate a problem in the preprocessing chain.
    """
    raise NotImplementedError


def _align_and_trim(
    charge_wfs: list[np.ndarray],
    current_wfs: list[np.ndarray],
    time_axes: list[np.ndarray],
    tp_aoe_max_values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Align waveforms to ``tp_aoe_max`` and trim to the common overlapping window.

    For each waveform, the time axis is shifted by subtracting its
    ``tp_aoe_max``, so that zero corresponds to the current maximum. The
    common overlapping window is then::

        t_min = max(t.min() for t in shifted_time_axes)
        t_max = min(t.max() for t in shifted_time_axes)

    and all waveforms are trimmed to this window.

    This function differs from ``psl.py:align_waveforms_to_peak`` in three
    ways:

    - Uses pre-computed ``tp_aoe_max_values`` rather than finding the argmax
      internally (the argmax was already computed on the MWA'd current in
      ``_recalculate_tp_aoe_max``; recomputing it here on the charge waveform
      would give a different result).
    - Works in physical time units (ns) rather than sample indices.
    - Determines the output window dynamically from the data rather than
      requiring a fixed output length.

    Parameters
    ----------
    charge_wfs : list[np.ndarray]
        List of normalised, upsampled charge waveforms.
    current_wfs : list[np.ndarray]
        List of MWA'd current waveforms. Must correspond element-wise to
        ``charge_wfs`` and share the same time axes.
    time_axes : list[np.ndarray]
        List of time axes in ns (one per waveform, not yet aligned).
    tp_aoe_max_values : np.ndarray
        1D array of recalculated ``tp_aoe_max`` values in ns, one per
        waveform, as returned by ``_recalculate_tp_aoe_max``.

    Returns
    -------
    common_times : np.ndarray
        1D array of shape ``(n_samples,)`` in ns. Zero = ``tp_aoe_max``.
    aligned_charge_wfs : np.ndarray
        2D array of shape ``(n_events, n_samples)``.
    aligned_current_wfs : np.ndarray
        2D array of shape ``(n_events, n_samples)``.
    """
    raise NotImplementedError


# ---------------------------------------------------------------------------
# Superpulse computation (one pixel)
# ---------------------------------------------------------------------------


def compute_superpulse(
    common_times: np.ndarray,
    charge_wfs: np.ndarray,
    current_wfs: np.ndarray,
    pixel: Pixel,
    detector: str,
    n_events_preliminary: int,
) -> Superpulse:
    """
    Compute the average charge and current waveforms for one pixel.

    Called twice per pixel: once on all preprocessed waveforms (preliminary
    superpulse), and once on the golden subset after ``apply_chi2_cut``
    (final superpulse).

    Parameters
    ----------
    common_times : np.ndarray
        Common time axis in ns, shape ``(n_samples,)``.
    charge_wfs : np.ndarray
        2D array of shape ``(n_events, n_samples)``.
    current_wfs : np.ndarray
        2D array of shape ``(n_events, n_samples)``.
    pixel : Pixel
        The pixel this superpulse belongs to.
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
    raise NotImplementedError


def compute_chi2_vs_superpulse(
    charge_wfs: np.ndarray,
    superpulse: Superpulse,
    baseline_region_mask: np.ndarray | None = None,
) -> np.ndarray:
    """
    Compute the reduced chi-squared of each charge waveform vs the superpulse.

    The noise sigma is estimated from the standard deviation of each waveform
    in the baseline region (``superpulse.time_axis < 0``, i.e. before
    ``tp_aoe_max``). A separate sigma is estimated for the superpulse from
    the same region and added in quadrature.

    Parameters
    ----------
    charge_wfs : np.ndarray
        2D array of shape ``(n_events, n_samples)``.
    superpulse : Superpulse
        Preliminary superpulse from ``compute_superpulse``.
    baseline_region_mask : np.ndarray or None, optional
        Boolean mask of shape ``(n_samples,)`` selecting the baseline region.
        If ``None``, derived from ``superpulse.time_axis < 0``.

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

    where ``n_dof = n_samples``. The use of the baseline region for sigma
    estimation is a practical approximation; a more principled approach would
    use an independent noise measurement. This is a known limitation.
    """
    raise NotImplementedError


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
        back to the original ``pixel_data_with_wfs``.
    """
    raise NotImplementedError


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def build_superpulse_for_pixel(
    pixel: Pixel,
    detector: str,
    tab_map: dict[str, int],
    evt_files: list[str],
    dsp_files: list[str],
    raw_files: list[str],
    dsp_fields: list[str],
    bl_window: tuple[float, float],
    n_target_wfs: int = 100,
    chi2_threshold: float = 3.0,
    upsample_factor: int = 16,
    upsample_mode: str = "l",
    mwa_length: int = 48,
    mwa_num: int = 3,
    mwa_type: int = 0,
) -> Superpulse:
    """
    Build the final superpulse for one pixel using an adaptive two-phase loop.

    This is the main entry point for superpulse construction. It manages the
    file loop and the two-phase adaptive strategy:

    **Phase 1 - accumulate:** read files one at a time, running the full
    dataset preparation pipeline (steps 1-7) on each. Accumulate the
    resulting ``pixel_data_with_wfs`` arrays. Continue until the accumulated
    event count reaches ``n_target_wfs``.

    **Phase 2 - build and check:** run the full superpulse construction
    pipeline (steps 8-12) on the accumulated data. If the number of golden
    waveforms after the chi2 cut is below ``n_target_wfs``, return to
    Phase 1 and continue reading from the next unread file, adding to the
    existing accumulation. Repeat until either:

    - ``n_golden >= n_target_wfs``: success, return the final superpulse.
    - All files are exhausted: log a warning and return the best superpulse
      available, or raise if no waveforms were found at all.

    The preprocessing (upsampling, MWA, alignment) is re-run on the full
    accumulated set each time Phase 2 is entered, not just on the new events.
    This ensures the common time window and alignment are globally consistent.

    Parameters
    ----------
    pixel : Pixel
        The energy-drift-time pixel to process.
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
    dsp_fields : list[str]
        DSP fields to attach to events, e.g.
        ``["tp_aoe_max", "cuspEmax", "tp_0_est"]``.
        ``tp_aoe_max`` must be included.
    bl_window : tuple[float, float]
        Baseline window ``(t_start, t_end)`` in ns, in the original
        non-aligned time axis.
    n_target_wfs : int, optional
        Target number of golden waveforms (after chi2 cut) for the final
        superpulse. Default 100.
    chi2_threshold : float, optional
        Reduced chi-squared threshold for the self-similarity cut. Default 3.0.
    upsample_factor : int, optional
        Upsampling factor. Default 16.
    upsample_mode : str, optional
        Interpolation mode for upsampling. Default ``'l'`` (linear).
    mwa_length : int, optional
        Moving window length in samples. Default 48.
        Must match ``psl.py:make_realistic_pulse_shape_lib``.
    mwa_num : int, optional
        Number of moving window passes. Default 3.
        Must match ``psl.py:make_realistic_pulse_shape_lib``.
    mwa_type : int, optional
        Moving window type. Default 0. Must match ``psl.py``.

    Returns
    -------
    Superpulse
        Final superpulse with ``n_events_final >= n_target_wfs`` if enough
        data was available, otherwise the best achievable superpulse from all
        available files.

    Raises
    ------
    ValueError
        If no waveforms are found in the pixel across all files.
    RuntimeError
        If ``evt_files``, ``dsp_files``, and ``raw_files`` have different
        lengths.

    Notes
    -----
    File triplets ``(evt_files[i], dsp_files[i], raw_files[i])`` must
    correspond to the same run segment. The caller is responsible for
    ensuring this ordering.

    Memory usage is proportional to the number of accumulated events times
    the waveform length times the upsampling factor. For typical pixel sizes
    (O(200) events, O(4000) samples at 16 ns, upsampled to O(64000) samples
    at 1 ns) this is of order a few hundred MB, which is acceptable.
    """
    raise NotImplementedError


# ---------------------------------------------------------------------------
# Output: all pixels, one detector
# ---------------------------------------------------------------------------


def write_superpulses_to_lh5(
    superpulses: dict[Pixel, Superpulse],
    output_path: str,
    detector: str,
) -> None:
    """
    Write all per-pixel superpulses for one detector to a single LH5 file.

    Iterates over ``superpulses``, calls ``Superpulse.to_lgdo()`` for each,
    and writes the result as a named group inside the file.

    The LH5 file structure is::

        {detector}/
            dt_{lo}_{hi}_ns/
                charge_wf             [Array, n_samples]
                current_wf            [Array, n_samples]
                time_axis             [Array, n_samples, units=ns]
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
    superpulses : dict[Pixel, Superpulse]
        Dictionary mapping each pixel to its final superpulse, as accumulated
        in the per-pixel processing loop.
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
    raise NotImplementedError