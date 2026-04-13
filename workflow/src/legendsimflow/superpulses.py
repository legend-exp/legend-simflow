"""
This module implements the dataset preparation and average ("superpulse") 
construction to characterize the pulse shape response of HPGe detectors.

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
            time_axis=t,
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
        time_axis: np.ndarray,
        slice: Slice,
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
        DSP fields to attach, e.g. ``["tp_aoe_max", "tp_0_est"]``.
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
    each event. For each channel in ``tab_map``, it identifies the corresponding
    TCM rows, reads only those rows from the DSP file via 
    ``lgdo.lh5.read(..., idx=sorted_idx)``, and then re-orders results to match 
    the original event ordering.

    ``lh5.read()`` requires a sorted ``idx`` array. The default behavior
    (``use_h5idx=False``) reads the full column and then indexes in memory. This 
    is faster than random HDF5 row access for most sizes, but uses more memory. 
    For very sparse reads, this trade-off should be re-evaluated.

    Raises
    ------
    ValueError
        If ``tp_aoe_max`` is not in ``fields`` (required for drift_time), or if 
        any rawid in the event channels is absent from ``tab_map``.
    """
    raise NotImplementedError


def _read_dsp_field_for_channel(
    channels: ak.Array,
    rows: ak.Array,
    dsp_file: str,
    field: str,
    tab_map: dict[str, int],
) -> ak.Array:
    """Read one DSP field for all channels in a set of events using TCM indices."""
    raise NotImplementedError


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
    raise NotImplementedError


# ---------------------------------------------------------------------------
# Waveform extraction (one slice)
# ---------------------------------------------------------------------------


def get_charge_and_current_wfs_for_slice(
    raw_file: Path | str,
    lh5_group: str,
    indices: list[int],
    dsp_config: Path | str,
    charge_output: str = "wf_blsub",
    current_output: str = "curr_av",
    align: str = "tp_aoe_max",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract aligned charge and current waveforms for a set of slice events,
    using the production DSP processing chain via ``WaveformBrowser``.

    A single ``WaveformBrowser`` instance is created for the file, with both
    ``charge_output`` and ``current_output`` requested as ``lines``. Events
    are iterated with ``find_entry`` and the processed waveforms are extracted
    from the browser's internal ``Line2D`` objects and stacked into arrays.

    The ``dsp_config`` and ``lh5_group`` required by this function can be
    obtained from ``legendsimflow.hpge_pars.lookup_currmod_fit_inputs``, which
    resolves both from the L200 data production directory given a run ID and
    detector name. The ``indices`` correspond to the raw-tier row indices for
    the slice events, which can be derived from the ``raw_wf_pairs`` it
    returns.

    Parameters
    ----------
    raw_file : Path or str
        Path to the raw-tier LH5 file containing the waveforms.
    lh5_group : str
        HDF5 group containing the waveform table, e.g. ``"ch1084803/raw"``.
    indices : list[int]
        Raw-tier row indices of the slice events to extract. These must all
        belong to ``raw_file``. Derived from the ``raw_wf_pairs`` returned by
        ``legendsimflow.hpge_pars.lookup_currmod_fit_inputs``.
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
    times : np.ndarray
        1D array of shape ``(n_samples,)`` in ns. Common time axis, aligned
        so that ``tp_aoe_max = 0``. Taken from the first valid event; all
        events share the same axis after alignment.
    charge_wfs : np.ndarray
        2D array of shape ``(n_valid_events, n_samples)``. Baseline-subtracted
        charge waveforms, normalised by ``cuspEmax`` so the plateau equals 1,
        aligned to ``tp_aoe_max``.
    current_wfs : np.ndarray
        2D array of shape ``(n_valid_events, n_samples)``. MWA'd current
        waveforms aligned to ``tp_aoe_max``.

    Notes
    -----
    A single ``WaveformBrowser`` instance is reused across all events in
    ``indices`` to avoid reconstructing the ``dspeed`` ``ProcessingChain``
    on every call. ``find_entry`` is called once per event in a loop;
    waveform data is extracted via ``browser.lines[output][0].get_xdata()``
    and ``.get_ydata()``.

    The normalization parameter ``cuspEmax`` is evaluated on the fly by the 
    DSP chain alongside the requested waveforms.
    Events for which the browser returns NaN-valued waveforms are dropped. 
    The returned arrays contain only valid events.

    Raises
    ------
    ValueError
        If no valid waveforms are found across all ``indices``.
    """
    # HACK: deferred import to avoid pint unit registry conflict at module level
    from dspeed.vis import WaveformBrowser  # noqa: PLC0415

    raise NotImplementedError

# ---------------------------------------------------------------------------
# Superpulse computation (one slice)
# ---------------------------------------------------------------------------


def compute_superpulse(
    common_times: np.ndarray,
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
    common_times : np.ndarray
        Common time axis in ns, shape ``(n_samples,)``.
    charge_wfs : np.ndarray
        2D array of shape ``(n_events, n_samples)``.
    current_wfs : np.ndarray
        2D array of shape ``(n_events, n_samples)``.
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
    raise NotImplementedError


def compute_chi2_vs_superpulse(
    charge_wfs: np.ndarray,
    superpulse: Superpulse,
    baseline_region_mask: np.ndarray | None = None,
) -> np.ndarray:
    """
    Compute the reduced chi-squared of each charge waveform vs the superpulse.

    The noise sigma is estimated from the standard deviation of each waveform
    in the baseline region. A separate sigma is estimated for the superpulse from
    the same region and added in quadrature.

    Parameters
    ----------
    charge_wfs : np.ndarray
        2D array of shape ``(n_events, n_samples)``.
    superpulse : Superpulse
        Preliminary superpulse from ``compute_superpulse``.
    baseline_region_mask : np.ndarray or None, optional
        Boolean mask of shape ``(n_samples,)`` selecting the baseline region.
        If ``None``, derived from ``superpulse.time_axis < (0 - drift_time)``. 

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
    use an independent noise measurement. 
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
        back to the original ``slice_data_with_wfs``.
    """
    raise NotImplementedError


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
    charge_output: str = "wf_blsub",
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
    raise NotImplementedError


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
    raise NotImplementedError