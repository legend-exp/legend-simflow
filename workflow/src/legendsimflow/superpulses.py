"""
Module to implement the dataset preparation and average ("superpulse") construction to characterize the pulse shape response of HPGe detectors.

This is an important step in tuning the pulse shape discrimination (PSD) simulations

Conventions
-----------
- All times are in nanoseconds and energies in keV unless otherwise stated.
- The ``tab_map`` convention follows the format:
  ``{detector_name: rawid}`` e.g. ``{"V03422A": 1084803}``
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path

import awkward as ak
import numpy as np
from legendmeta import LegendMetadata
from lgdo import Array, Scalar, Struct, lh5
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from legendsimflow import metadata as mutils
from legendsimflow import utils

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
    ):
        """
        Serialise to an ``Struct`` ready for writing to LH5.

        Returns
        -------
        Struct
            Struct with the following fields:

            - ``charge_wf``            : ``Array``, shape ``(n_charge_samples,)``
            - ``current_wf``           : ``Array``, shape ``(n_current_samples,)``
            - ``charge_time_axis``     : ``Array``, shape ``(n_charge_samples,)``, attrs ``{"units": "ns"}``
            - ``current_time_axis``    : ``Array``, shape ``(n_current_samples,)``, attrs ``{"units": "ns"}``
            - ``dt_center``            : ``Scalar``, drift time center [ns]
            - ``dt_lo``                : ``Scalar``, drift time lower bound [ns]
            - ``dt_hi``                : ``Scalar``, drift time upper bound [ns]
            - ``e_lo``                 : ``Scalar``, energy lower bound [keV]
            - ``e_hi``                 : ``Scalar``, energy upper bound [keV]
            - ``detector``             : ``Scalar``, detector name string
            - ``n_events_preliminary`` : ``Scalar``
            - ``n_events_final``       : ``Scalar``
        """
        return Struct(
            {
                "charge_wf": Array(self.charge_wf),
                "current_wf": Array(self.current_wf),
                "charge_time_axis": Array(self.charge_time_axis, attrs={"units": "ns"}),
                "current_time_axis": Array(
                    self.current_time_axis, attrs={"units": "ns"}
                ),
                "dt_center": Scalar(
                    self.slice.drift_time_center, attrs={"units": "ns"}
                ),
                "dt_lo": Scalar(self.slice.drift_time_range[0], attrs={"units": "ns"}),
                "dt_hi": Scalar(self.slice.drift_time_range[1], attrs={"units": "ns"}),
                "e_lo": Scalar(self.slice.energy_range[0], attrs={"units": "keV"}),
                "e_hi": Scalar(self.slice.energy_range[1], attrs={"units": "keV"}),
                "detector": Scalar(self.detector),
                "n_events_preliminary": Scalar(self.n_events_preliminary),
                "n_events_final": Scalar(self.n_events_final),
            }
        )


def lookup_superpulse_inputs(
    l200data: str | Path,
    metadata: LegendMetadata,
    runid: str,
    hpge: str,
    max_files: int | None = None,
) -> tuple[list[Path], list[Path], list[Path], Path, dict[str, int]]:
    """Look up all inputs needed to build superpulses for one detector and run.

    Parameters
    ----------
    l200data
        Path to the L200 data production directory.
    metadata
        The metadata instance.
    runid
        LEGEND run identifier, e.g. ``"l200-p16-r008-ssc"``.
    hpge
        Detector name, e.g. ``"V03422A"``.
    max_files
        Limit the number of files per tier. Default: all.

    Returns
    -------
    raw_files, dsp_files, evt_files
        Sorted lists of LH5 file paths for each tier.
    dsp_cfg_file
        Path to the DSP configuration file.
    tab_map
        Mapping ``{detector_name: rawid}``.
    """
    if isinstance(l200data, str):
        l200data = Path(l200data)

    _, period_int, run_int, data_type = mutils.parse_runid(runid)
    period = f"p{period_int:02d}"
    run = f"r{run_int:03d}"
    df_cfg = utils.lookup_dataflow_config(l200data).paths

    raw_files = sorted((df_cfg.tier_raw / data_type / period / run).glob("*.lh5"))[
        :max_files
    ]
    dsp_files = sorted((df_cfg.tier_dsp / data_type / period / run).glob("*.lh5"))[
        :max_files
    ]
    evt_files = sorted((df_cfg.tier_evt / data_type / period / run).glob("*.lh5"))[
        :max_files
    ]

    if not evt_files:
        msg = f"no EVT files found for {runid}"
        raise FileNotFoundError(msg)

    dsp_cfg_regex = r"l200-*-r%-T%-ICPC-dsp_proc_chain.*"
    dsp_cfg_files = list(
        (l200data / "inputs/dataprod/config/tier_dsp").glob(dsp_cfg_regex)
    )
    if dsp_cfg_files == []:
        dsp_cfg_files = list(
            (l200data / "inputs/dataprod/config/tier/dsp").glob(dsp_cfg_regex)
        )
    if len(dsp_cfg_files) != 1:
        msg = f"could not find a suitable dsp config file in {l200data} (or multiple found)"
        raise RuntimeError(msg)

    tstamp = mutils.runinfo(metadata, runid).start_key
    chmap = metadata.channelmap(tstamp)
    tab_map = {hpge: chmap[hpge]["daq"]["rawid"]}

    return raw_files, dsp_files, evt_files, dsp_cfg_files[0], tab_map


def read_and_select_evt_data(
    evt_files: str | list[str],
    aoe_low_threshold: float = -3.0,
    aoe_high_threshold: float = 3.0,
) -> ak.Array:
    """Read EVT-tier LH5 file(s) and apply quality and PSD cuts.

    Only the fields needed for event selection are read.

    Parameters
    ----------
    evt_files
        Path to one or more EVT-tier LH5 files.
    aoe_low_threshold
        Lower bound on ``geds.psd.low_aoe.value``. Default ``-3.0``.
    aoe_high_threshold
        Upper bound on ``geds.psd.high_aoe.value``. Default ``3.0``.

    Returns
    -------
        Filtered event data passing all quality and PSD cuts.
    """
    evt_data = lh5.read_as(
        "evt",
        evt_files,
        library="ak",
        field_mask=[
            "geds",
            "coincident",
            "trigger",
            "spms/energy_sum",
            "spms/first_t0",
        ],
    )

    mask = (
        ak.all(evt_data.geds.quality.is_good_channel, axis=-1)
        & ~evt_data.trigger.is_forced
        & ~evt_data.coincident.puls
        & ~evt_data.coincident.muon
        & ~evt_data.coincident.muon_offline
        & evt_data.geds.quality.is_bb_like
        & (evt_data.geds.multiplicity == 1)
        & (evt_data.spms.energy_sum > 10)
    )
    evt_data = evt_data[mask]

    psd_mask = ak.all(
        evt_data.geds.psd.low_aoe.value > aoe_low_threshold, axis=-1
    ) & ak.all(evt_data.geds.psd.high_aoe.value < aoe_high_threshold, axis=-1)

    return evt_data[psd_mask]


def select_detector_events(
    evt_data: ak.Array,
    detector: str,
) -> ak.Array:
    """
    Filter event data to only events from a single detector.

    Called after ``read_and_select_evt_data`` and before
    ``add_dsp_pars_to_evt``, so that the DSP read operates on the smallest
    possible array.

    Parameters
    ----------
    evt_data
        Event data after quality and PSD cuts, shape ``(n_events,)``.
    detector
        Name of the target detector, e.g. ``"V03422A"``.

    Returns
    -------
        Filtered event data containing only events where
        ``geds.detector_name == detector``, shape ``(n_det_events,)``.
    """
    mask = ak.any(evt_data.geds.detector_name == detector, axis=-1)
    det_evt_data = evt_data[mask]

    log.debug(
        "detector %s: %d -> %d events",
        detector,
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
    det_evt_data
        Single-detector event data as returned by ``select_detector_events``,
        shape ``(n_det_events,)``.
    dsp_file
        DSP-tier LH5 file path, corresponding to the same run segment as the
        EVT file used to produce ``det_evt_data``. Must be a single file
        because ``hit_idx`` values are row indices into this specific file.
    tab_map
        Mapping from detector name to rawid. Only the entry for ``detector``
        is used.
    fields
        DSP fields to attach, e.g. ``["tp_aoe_max", "tp_0_est"]``.
        ``tp_aoe_max`` must be included for ``drift_time`` to be computed.

    Returns
    -------
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
    det_evt_data
        Single-detector event data with DSP fields attached, as returned by
        ``add_dsp_pars_to_evt``, shape ``(n_det_events,)``.
    slice
        The energy-drift-time slice to select.

    Returns
    -------
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


def get_charge_and_current_wfs_for_slice(
    raw_file: Path | str,
    lh5_group: str,
    indices: list[int],
    dsp_config: Path | str,
    charge_output: str = "wf_pz_win",
    current_output: str = "curr_av",
    align: str = "tp_aoe_max",
) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray | None, np.ndarray | None
]:
    """
    Extract aligned charge and current waveforms for a set of slice events.

    Uses the DSP processing chain via ``WaveformBrowser``. After alignment each
    event has the same sampling rate and number of samples but a different
    x-offset. This function collects all events, finds the overlapping time
    region, and trims every waveform to that common window.

    Parameters
    ----------
    raw_file
        Path to the raw-tier LH5 file.
    lh5_group
        HDF5 group containing the waveform table, e.g. ``"ch1084803/raw"``.
    indices
        Raw-tier row indices of the slice events.
    dsp_config
        Path to the production DSP configuration JSON file.
    charge_output
        DSP output name for the charge waveform.
    current_output
        DSP output name for the current waveform.
    align
        DSP parameter used to align waveforms on the time axis.

    Returns
    -------
    charge_times
        Common time axis for charge waveforms, shape ``(n_common_charge,)``.
    current_times
        Common time axis for current waveforms, shape ``(n_common_current,)``.
    charge_wfs
        Trimmed charge waveforms, shape ``(n_valid, n_common_charge)``.
    current_wfs
        Trimmed current waveforms, shape ``(n_valid, n_common_current)``.
    bl_std
        Baseline std (ADC) per event, shape ``(n_valid,)``.
    cuspEmax
        Energy estimator (ADC) per event, shape ``(n_valid,)``.

    Raises
    ------
    ValueError
        If ``indices`` is empty or no valid waveforms are found.
    """
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
        n_drawn=len(indices),
    )

    browser.find_next(n_wfs=len(indices), append=False)

    charge_lines = browser.lines.get(charge_output, [])
    current_lines = browser.lines.get(current_output, [])

    if not charge_lines or not current_lines:
        msg = f"no waveforms returned from {raw_file}"
        raise ValueError(msg)

    bl_std_vals = browser.legend_vals.get("bl_std", [])
    cuspEmax_vals = browser.legend_vals.get("cuspEmax", [])

    # First pass: collect valid events with their x- and y-data
    valid = []
    for i, (cl, il) in enumerate(zip(charge_lines, current_lines, strict=True)):
        charge_y = cl.get_ydata()
        current_y = il.get_ydata()

        if np.any(np.isnan(charge_y)) or np.any(np.isnan(current_y)):
            log.debug("event %d: NaN in waveform, skipping", i)
            continue

        entry = {
            "charge_x": cl.get_xdata(),
            "charge_y": charge_y,
            "current_x": il.get_xdata(),
            "current_y": current_y,
        }
        if bl_std_vals and cuspEmax_vals:
            entry["bl_std"] = float(bl_std_vals[i].magnitude)
            entry["cuspEmax"] = float(cuspEmax_vals[i].magnitude)

        valid.append(entry)

    if not valid:
        msg = f"no valid waveforms found in {raw_file} for {len(indices)} indices"
        raise ValueError(msg)

    # Find the common time window (max of left edges, min of right edges)
    charge_t_min = max(e["charge_x"][0] for e in valid)
    charge_t_max = min(e["charge_x"][-1] for e in valid)
    current_t_min = max(e["current_x"][0] for e in valid)
    current_t_max = min(e["current_x"][-1] for e in valid)

    # Trim each event to the common window
    charge_list = []
    current_list = []
    bl_std_list = []
    cuspEmax_list = []

    for e in valid:
        c_mask = (e["charge_x"] >= charge_t_min) & (e["charge_x"] <= charge_t_max)
        i_mask = (e["current_x"] >= current_t_min) & (e["current_x"] <= current_t_max)

        charge_list.append(e["charge_y"][c_mask])
        current_list.append(e["current_y"][i_mask])

        if "bl_std" in e:
            bl_std_list.append(e["bl_std"])
            cuspEmax_list.append(e["cuspEmax"])

    # Common time axes from the first event's trimmed region
    e0 = valid[0]
    charge_times = e0["charge_x"][
        (e0["charge_x"] >= charge_t_min) & (e0["charge_x"] <= charge_t_max)
    ]
    current_times = e0["current_x"][
        (e0["current_x"] >= current_t_min) & (e0["current_x"] <= current_t_max)
    ]

    log.debug(
        "extracted %d/%d valid waveforms from %s",
        len(valid),
        len(indices),
        Path(raw_file).name,
    )

    return (
        charge_times,
        current_times,
        np.array(charge_list),
        np.array(current_list),
        np.array(bl_std_list) if bl_std_list else None,
        np.array(cuspEmax_list) if cuspEmax_list else None,
    )


def compute_superpulse(
    charge_times: np.ndarray,
    current_times: np.ndarray,
    charge_wfs: np.ndarray,
    current_wfs: np.ndarray,
    slice: Slice,
    detector: str,
    n_events: int,
) -> Superpulse:
    """
    Compute the average charge and current waveforms for one slice.

    Called twice per slice: once on all preprocessed waveforms (preliminary
    superpulse), and once on the golden subset after ``apply_chi2_cut``
    (final superpulse).

    Parameters
    ----------
    charge_times
        Time axis for the charge waveforms in ns, shape ``(n_charge_samples,)``.
    current_times
        Time axis for the current waveforms in ns, shape ``(n_current_samples,)``.
    charge_wfs
        2D array of shape ``(n_events, n_charge_samples)``.
    current_wfs
        2D array of shape ``(n_events, n_current_samples)``.
    slice
        The slice this superpulse belongs to.
    detector
        Detector name.
    n_events
        Total number of events.

    Returns
    -------
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
        n_events_preliminary=n_events,
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
    charge_wfs
        2D array of shape ``(n_events, n_samples)``.
    superpulse : Superpulse
        Preliminary superpulse from ``compute_superpulse``.
    baseline_region_mask
        Boolean mask of shape ``(n_samples,)`` selecting the baseline region.
        If ``None``, derived from
        ``superpulse.charge_time_axis < -superpulse.slice.drift_time_range[1]``.
        Only used in fallback mode.
    bl_std
        Per-event baseline standard deviation in ADC units, shape
        ``(n_events,)``. From ``get_charge_and_current_wfs_for_slice``.
    cuspEmax
        Per-event energy estimator in ADC units, shape ``(n_events,)``.
        From ``get_charge_and_current_wfs_for_slice``.

    Returns
    -------
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
        sigma_wfs = np.nanstd(
            charge_wfs[:, baseline_region_mask], axis=1
        )  # (n_events,)
        sigma_sp = np.nanstd(sp_wf[baseline_region_mask])
        if not np.all(np.isfinite(sigma_wfs)) or not np.isfinite(sigma_sp):
            msg = (
                "baseline noise estimation produced non-finite values; "
                "ensure the baseline_region_mask includes finite samples, "
                "or provide bl_std and cuspEmax from the DSP chain"
            )
            raise ValueError(msg)

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
    charge_wfs
        2D array of shape ``(n_events, n_samples)``.
    current_wfs
        2D array of shape ``(n_events, n_samples)``.
    chi2_values
        1D array of reduced chi-squared values from
        ``compute_chi2_vs_superpulse``.
    threshold
        Reduced chi-squared cut value. Default 3.0.

    Returns
    -------
    golden_charge_wfs
        Shape ``(n_golden, n_samples)``.
    golden_current_wfs
        Shape ``(n_golden, n_samples)``.
    golden_indices
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


def trim_and_stack(
    times_list: list[np.ndarray], wfs_list: list[np.ndarray]
) -> tuple[np.ndarray, np.ndarray]:
    """Find the common time window across multiple arrays, trim, and stack them."""
    tmin = max(t[0] for t in times_list)
    tmax = min(t[-1] for t in times_list)

    masks = [(t >= tmin) & (t <= tmax) for t in times_list]

    # Enforce minimum length to handle ±1 sample sub-grid alignment shifts
    nmin = min(m.sum() for m in masks)

    stacked_wfs = np.concatenate(
        [w[:, np.where(m)[0][:nmin]] for m, w in zip(masks, wfs_list, strict=True)],
        axis=0,
    )

    common_times = times_list[0][np.where(masks[0])[0][:nmin]]

    return common_times, stacked_wfs


def accumulate_wfs_for_slice(
    slice: Slice,
    detector: str,
    tab_map: dict[str, int],
    evt_files: list[str],
    dsp_files: list[str],
    raw_files: list[str],
    dsp_config: Path | str,
    dsp_fields: list[str],
    start_file_idx: int = 0,
    n_target: int = 150,
    charge_output: str = "wf_pz_win",
    current_output: str = "curr_av",
    align: str = "tp_aoe_max",
) -> tuple[
    np.ndarray | None,
    np.ndarray | None,
    np.ndarray | None,
    np.ndarray | None,
    np.ndarray | None,
    np.ndarray | None,
    int,
]:
    """Accumulate aligned waveforms for one slice, dynamically trimming to a global time window."""
    if not (len(evt_files) == len(dsp_files) == len(raw_files)):
        msg = "File lists must have equal lengths."
        raise RuntimeError(msg)
    if detector not in tab_map:
        msg = f"Unknown detector {detector!r}"
        raise KeyError(msg)

    lh5_group = f"ch{tab_map[detector]}/raw"
    c_times, i_times, c_wfs, i_wfs, bl_stds, cusps = [], [], [], [], [], []

    next_idx = start_file_idx
    for file_idx in range(start_file_idx, len(evt_files)):
        try:
            evt_data = read_and_select_evt_data(evt_files[file_idx])
            det_data = select_detector_events(evt_data, detector)
            if len(det_data) == 0:
                continue

            det_data = add_dsp_pars_to_evt(
                det_data, dsp_files[file_idx], tab_map, dsp_fields
            )
            slice_data = select_data_in_slice(det_data, slice)
            if len(slice_data) == 0:
                continue

            ct, it, cw, iw, bl, ce = get_charge_and_current_wfs_for_slice(
                raw_file=raw_files[file_idx],
                lh5_group=lh5_group,
                indices=ak.flatten(slice_data.geds.hit_idx).to_list(),
                dsp_config=dsp_config,
                charge_output=charge_output,
                current_output=current_output,
                align=align,
            )

            c_times.append(ct)
            i_times.append(it)
            c_wfs.append(cw)
            i_wfs.append(iw)
            if bl is not None:
                bl_stds.append(bl)
            if ce is not None:
                cusps.append(ce)

        except ValueError:
            log.debug("File %d: No valid waveforms extracted.", file_idx)
            continue
        finally:
            next_idx = file_idx + 1

        # Check target quota
        if sum(len(c) for c in c_wfs) >= n_target:
            break

    # Return Nones if no valid waveforms were found at all
    if not c_wfs:
        return None, None, None, None, None, None, next_idx

    c_common, c_stack = trim_and_stack(c_times, c_wfs)
    i_common, i_stack = trim_and_stack(i_times, i_wfs)

    bl_stack = np.concatenate(bl_stds) if bl_stds else None
    ce_stack = np.concatenate(cusps) if cusps else None

    log.info("%s | %s | Accumulated %d waveforms.", detector, slice, len(c_stack))

    return c_common, i_common, c_stack, i_stack, bl_stack, ce_stack, next_idx


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
    superpulses
        Dictionary mapping each slice to its final superpulse, as accumulated
        in the per-slice processing loop.
    output_path
        Path to the output LH5 file, e.g.
        ``"output/V03422A_superpulses.lh5"``.
    detector
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


def read_superpulses_from_lh5(
    path: str,
    detector: str,
) -> dict[Slice, Superpulse]:
    """Read superpulses written by :func:`write_superpulses_to_lh5`.

    Parameters
    ----------
    path
        Path to the LH5 file.
    detector
        Detector name (top-level group), e.g. ``"V03422A"``.

    Returns
    -------
    dict[Slice, Superpulse]
    """
    superpulses: dict[Slice, Superpulse] = {}

    # List subgroups under {detector}/
    det_struct = lh5.read(detector, path)

    for key in det_struct:
        # Keys are like "dt_1100_1300_ns"
        match = re.match(r"dt_(\d+)_(\d+)_ns", key)
        if match is None:
            log.debug("skipping unrecognised key %s/%s", detector, key)
            continue

        group = lh5.read(f"{detector}/{key}", path)

        sl = Slice(
            energy_range=(float(group["e_lo"].value), float(group["e_hi"].value)),
            drift_time_range=(float(group["dt_lo"].value), float(group["dt_hi"].value)),
        )

        sp = Superpulse(
            charge_wf=group["charge_wf"].view_as("np"),
            current_wf=group["current_wf"].view_as("np"),
            charge_time_axis=group["charge_time_axis"].view_as("np"),
            current_time_axis=group["current_time_axis"].view_as("np"),
            slice=sl,
            detector=group["detector"].value,
            n_events_preliminary=int(group["n_events_preliminary"].value),
            n_events_final=int(group["n_events_final"].value),
        )

        superpulses[sl] = sp

    log.info("read %d slices for %s from %s", len(superpulses), detector, path)
    return superpulses


def plot_wfs_and_superpulse(
    charge_times: np.ndarray,
    current_times: np.ndarray,
    golden_charge_wfs: np.ndarray,
    golden_current_wfs: np.ndarray,
    superpulse: Superpulse,
):
    """
    Plot golden charge and current waveforms with the final superpulse.

    Two-panel figure: top = charge waveforms + superpulse,
    bottom = current waveforms + superpulse.
    Each panel uses its own time axis, but both are clipped to the
    intersection of the two time ranges for consistent visualization.

    Parameters
    ----------
    charge_times, current_times
        Time axes (may have different ranges).
    golden_charge_wfs, golden_current_wfs
        Golden waveform arrays after chi2 cut, shape ``(n_events, n_samples)``.
    superpulse
        The final superpulse to overlay.

    Returns
    -------
    fig : matplotlib.figure.Figure
    (ax_charge, ax_current) : tuple of Axes
    """
    sl = superpulse.slice

    # Shared x range: intersection of the two time axes
    x_min = max(charge_times[0], current_times[0])
    x_max = min(charge_times[-1], current_times[-1])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    fig.suptitle(
        f"{superpulse.detector} - {sl}  ({superpulse.n_events_final} events)",
        fontsize=14,
    )

    # Charge
    ax1.plot(
        charge_times, golden_charge_wfs.T, color="deepskyblue", alpha=0.3, linewidth=1
    )
    ax1.plot(
        [],
        [],
        color="deepskyblue",
        alpha=0.8,
        linewidth=1,
        label=f"Golden wfs ({golden_charge_wfs.shape[0]})",
    )
    ax1.plot(
        superpulse.charge_time_axis,
        superpulse.charge_wf,
        color="black",
        linewidth=2,
        label="superpulse",
    )
    ax1.set_ylabel("ADC / cuspEmax")
    ax1.set_xlim(x_min, x_max)
    ax1.legend()
    ax1.grid(alpha=0.3, linestyle="--")

    # Current
    ax2.plot(
        current_times, golden_current_wfs.T, color="deepskyblue", alpha=0.3, linewidth=1
    )
    ax2.plot(
        [],
        [],
        color="deepskyblue",
        alpha=0.8,
        linewidth=1,
        label=f"Golden wfs ({golden_current_wfs.shape[0]})",
    )
    ax2.plot(
        superpulse.current_time_axis,
        superpulse.current_wf,
        color="black",
        linewidth=2,
        label="superpulse",
    )
    ax2.set_xlabel("Time [ns]")
    ax2.set_ylabel("d(ADC/cuspEmax)/dt")
    ax2.set_xlim(x_min, x_max)
    ax2.legend()
    ax2.grid(alpha=0.3, linestyle="--")

    fig.tight_layout()
    return fig, (ax1, ax2)


def plot_chi2_cut(
    chi2_values: np.ndarray,
    chi2_threshold: float,
    times: np.ndarray,
    wfs: np.ndarray,
    final_superpulse: Superpulse,
    curve: str = "charge",
) -> tuple[Figure, tuple]:
    """
    Two-panel validation plot for the chi2 self-similarity cut.

    Left: histogram of reduced chi2 with threshold line.
    Right: all waveforms color-coded by cut status (blue=pass, red=fail) with the final superpulse overlaid.

    Parameters
    ----------
    chi2_values
        Reduced chi2 per waveform, shape ``(n_events,)``.
    chi2_threshold
        Cut threshold.
    times
        Time axis, shape ``(n_samples,)``.
    wfs
        All waveforms before cut, shape ``(n_events, n_samples)``.
    final_superpulse : Superpulse
        Superpulse built from golden waveforms (after cut).
    curve
        ``"charge"`` or ``"current"``. Default ``"charge"``.

    Returns
    -------
    fig : matplotlib.figure.Figure
    (ax_hist, ax_wfs) : tuple of Axes
    """
    if curve == "charge":
        sp_times = final_superpulse.charge_time_axis
        sp_wf = final_superpulse.charge_wf
        ylabel = "ADC / cuspEmax"
    else:
        sp_times = final_superpulse.current_time_axis
        sp_wf = final_superpulse.current_wf
        ylabel = "d(ADC/cuspEmax)/dt"

    sl = final_superpulse.slice
    mask_pass = chi2_values < chi2_threshold
    idx_fail = np.where(~mask_pass)[0]
    idx_pass = np.where(mask_pass)[0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f"{final_superpulse.detector} - {sl}", fontsize=14)

    # Left: histogram
    bw = chi2_threshold / 8
    bins = np.arange(0, chi2_values.max() + bw, bw)
    ax1.hist(
        chi2_values,
        bins=bins,
        alpha=0.7,
        color="grey",
        histtype="step",
        label="All events",
    )
    ax1.hist(
        chi2_values[mask_pass],
        bins=bins,
        alpha=0.9,
        color="deepskyblue",
        histtype="stepfilled",
        label=rf"$\chi^2_r$ < {chi2_threshold}",
    )
    ax1.axvline(chi2_threshold, color="red", linestyle="--", linewidth=1)
    ax1.set_xlabel(r"Reduced $\chi^2$")
    ax1.set_ylabel("Counts")
    ax1.legend()

    # Right: waveforms color-coded by pass/fail
    if len(idx_fail) > 0:
        ax2.plot(times, wfs[idx_fail].T, color="firebrick", alpha=0.8, linewidth=1)
        ax2.plot(
            [],
            [],
            color="firebrick",
            alpha=0.8,
            linewidth=1,
            label=f"Fail ({len(idx_fail)})",
        )
    if len(idx_pass) > 0:
        ax2.plot(times, wfs[idx_pass].T, color="deepskyblue", alpha=0.8, linewidth=1)
        ax2.plot(
            [],
            [],
            color="deepskyblue",
            alpha=0.8,
            linewidth=1,
            label=f"Pass ({len(idx_pass)})",
        )

    ax2.plot(sp_times, sp_wf, color="black", linewidth=2, label="Final superpulse")

    ax2.set_xlabel("Time [ns]")
    ax2.set_ylabel(ylabel)
    ax2.set_xlim(-1000, 500)
    ax2.grid(alpha=0.3, linestyle="--")
    ax2.legend(loc="upper left")

    fig.tight_layout()
    return fig, (ax1, ax2)


def plot_superpulses(
    lh5_file: str,
    detector: str,
    curve: str = "charge",
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
) -> tuple:
    """
    Plot superpulses from all drift-time slices, color-coded by drift time.

    Parameters
    ----------
    lh5_file
        Path to the LH5 file produced by :func:`write_superpulses_to_lh5`.
    detector
        Detector name (top-level group in the file).
    curve
        ``"charge"`` or ``"current"``. Default ``"charge"``.
    xlim
        x-axis limits in ns. Default ``(-5000, 5000)``.
    ylim
        y-axis limits. If ``None``, matplotlib auto-scales.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax  : matplotlib.axes.Axes
    """
    import matplotlib.colors as mcolors  # noqa: PLC0415
    from matplotlib import cm  # noqa: PLC0415

    if curve not in ("charge", "current"):
        msg = f"curve must be 'charge' or 'current', got {curve!r}"
        raise ValueError(msg)

    groups = sorted(
        lh5.ls(lh5_file, f"{detector}/"),
        key=lambda g: int(g.split("_")[1]),
    )

    norm = mcolors.Normalize(vmin=0, vmax=2500)
    viridis = cm.get_cmap("viridis")
    ylabel = "ADC / cuspEmax" if curve == "charge" else "d(ADC/cuspEmax)/dt"

    fig, ax = plt.subplots(figsize=(12, 6))

    e_lo, e_hi = None, None
    for group in groups:
        struct = lh5.read(group, lh5_file)

        if e_lo is None:
            e_lo = struct["e_lo"].value
            e_hi = struct["e_hi"].value

        dt_center = struct["dt_center"].value
        dt_lo = struct["dt_lo"].value
        dt_hi = struct["dt_hi"].value
        n_events = struct["n_events_final"].value
        times = struct[f"{curve}_time_axis"].nda
        wf = struct[f"{curve}_wf"].nda

        color = mcolors.to_hex(viridis(norm(dt_center)))
        ax.plot(
            times,
            wf,
            color=color,
            alpha=0.8,
            linewidth=2,
            label=f"dt=[{dt_lo:.0f}, {dt_hi:.0f}] ns  ({n_events} ev)",
        )

    ax.set_title(f"{detector} - {curve} superpulses  |  E = [{e_lo}, {e_hi}] keV")
    ax.set_xlabel("Time [ns]")
    ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.grid(visible=True, which="both", linestyle="--", alpha=0.5)

    sm = plt.cm.ScalarMappable(cmap=viridis, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label("Drift Time [ns]", rotation=270, labelpad=20)
    fig.tight_layout()

    return fig, ax
