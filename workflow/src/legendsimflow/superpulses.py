"""
Module to implement the dataset preparation and average ("superpulse") construction to characterize the pulse shape response of HPGe detectors.

This is an important step in tuning the pulse shape discrimination (PSD) simulations
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path

import awkward as ak
import dbetto
import numpy as np
from dbetto import AttrsDict
from legendmeta import LegendMetadata
from lgdo import Array, Scalar, Struct, lh5
from matplotlib import colormaps
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from legendsimflow import metadata as mutils
from legendsimflow import utils

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class Slice:
    """
    Defines a 2D slice in the energy-drift-time space.

    Parameters
    ----------
    energy_range
        Lower and upper bounds of the energy slice, in keV.
        Example: ``(1500.0, 2000.0)``
    drift_time_range
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
        """Construct the superpulse.

        Parameters
        ----------
        charge_wf
            Average normalised charge waveform, shape ``(n_charge_samples,)``.
            Amplitude is dimensionless (ADC / cuspEmax, normalised to 1 at
            plateau).
        current_wf
            Average current waveform, shape ``(n_current_samples,)``.
            Units: ``(ADC / cuspEmax) / ns * dt_data``, matching the convention
            in ``psl.py``.
        charge_time_axis
            Time axis for the charge waveform in ns, shape ``(n_charge_samples,)``,
            aligned so that ``tp_aoe_max = 0``.
        current_time_axis
            Time axis for the current waveform in ns, shape ``(n_current_samples,)``,
            aligned so that ``tp_aoe_max = 0``.
        slice
            The energy-drift-time slice this superpulse represents.
        detector
            Detector name, e.g. ``"V03422A"``.
        n_events_preliminary
            Number of waveforms used to build the preliminary superpulse (before
            the chi2 cut).
        n_events_final
            Number of waveforms surviving the chi2 cut, used to build this
            superpulse.
        """
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
            - ``drift_time_center``    : ``Scalar``, drift time center [ns]
            - ``drift_time_lo``        : ``Scalar``, drift time lower bound [ns]
            - ``drift_time_hi``        : ``Scalar``, drift time upper bound [ns]
            - ``energy_lo``            : ``Scalar``, energy lower bound [keV]
            - ``energy_hi``            : ``Scalar``, energy upper bound [keV]
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
                "drift_time_center": Scalar(
                    self.slice.drift_time_center, attrs={"units": "ns"}
                ),
                "drift_time_lo": Scalar(
                    self.slice.drift_time_range[0], attrs={"units": "ns"}
                ),
                "drift_time_hi": Scalar(
                    self.slice.drift_time_range[1], attrs={"units": "ns"}
                ),
                "energy_lo": Scalar(self.slice.energy_range[0], attrs={"units": "keV"}),
                "energy_hi": Scalar(self.slice.energy_range[1], attrs={"units": "keV"}),
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
    *,
    evt_tier_name: str = "evt",
) -> tuple[list[Path], list[Path], Path, dict[str, int], str]:
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
    evt_tier_name
        Name of the evt tier to look for, e.g. "evt" or "pet". Default: "evt".

    Returns
    -------
    raw_files, evt_files
        Sorted lists of LH5 file paths for the raw and evt tiers.
    dsp_cfg_file
        Path to the DSP configuration file.
    tab_map
        Mapping ``{detector_name: rawid}``.
    data_runid
        Run whose data was looked up: the reference calibration run for a
        physics ``runid``, otherwise ``runid`` itself.
    """
    if isinstance(l200data, str):
        l200data = Path(l200data)

    _, period_int, run_int, data_type = mutils.parse_runid(runid)

    # for physics runs we build superpulses from the reference calibration run:
    # its raw/evt data and channel map are what we actually read below
    data_runid = runid
    if data_type == "phy":
        data_runid = mutils.reference_cal_run(metadata, runid)
        _, period_int, run_int, data_type = mutils.parse_runid(data_runid)

        msg = f"inferred reference calibration run: {data_runid}"
        log.debug(msg)

    period = f"p{period_int:02d}"
    run = f"r{run_int:03d}"
    df_cfg = utils.lookup_dataflow_config(l200data).paths

    raw_files = sorted((df_cfg.tier_raw / data_type / period / run).glob("*.lh5"))[
        :max_files
    ]
    evt_files = sorted(
        (df_cfg[f"tier_{evt_tier_name}"] / data_type / period / run).glob("*.lh5")
    )[:max_files]

    if not evt_files:
        msg = f"no evt tier files found for {data_runid}."
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

    tstamp = mutils.runinfo(metadata, data_runid).start_key
    chmap = metadata.channelmap(tstamp, skip_version_check=True)
    tab_map = {hpge: chmap[hpge]["daq"]["rawid"]}

    return raw_files, evt_files, dsp_cfg_files[0], tab_map, data_runid


def _read_and_sel_evts(
    evt_files: str | list[str],
    detector: str,
    t0_field: str | None = None,
    aoe_low_threshold: float = -3.0,
    aoe_high_threshold: float = 3.0,
) -> ak.Array:
    """Read evt data and perform basic selections."""
    field_mask = ["geds", "coincident", "trigger"]
    if t0_field is not None:
        field_mask.append(t0_field)

    evt_data = lh5.read_as(
        "evt",
        evt_files,
        library="ak",
        field_mask=field_mask,
    )

    mask = (
        ak.all(evt_data.geds.quality.is_good_channel, axis=-1)
        & (~evt_data.trigger.is_forced)
        & (~evt_data.coincident.puls)
        & (~evt_data.coincident.muon)
        & (~evt_data.coincident.muon_offline)
        & evt_data.geds.quality.is_bb_like
        & (evt_data.geds.multiplicity == 1)
        & (ak.all(evt_data.geds.detector_name == detector, axis=-1))
    )

    # `t0_field` names the event start time used for the drift time; it is set
    # in the metadata settings. Drop events where it is NaN: for `spms/event_t0`
    # those are exactly the low-p.e. events. `t0_field` is None only when the
    # drift time is taken directly from `end_time_field`, in which case there is
    # no t0 to cut.
    if t0_field is not None:
        mask = mask & ~np.isnan(_get_nested_field(evt_data, t0_field))

    evt_data = evt_data[mask]

    psd_mask = ak.all(
        evt_data.geds.psd.low_aoe.value > aoe_low_threshold, axis=-1
    ) & ak.all(evt_data.geds.psd.high_aoe.value < aoe_high_threshold, axis=-1)

    return evt_data[psd_mask]


def _get_nested_field(data: ak.Array, field: str) -> ak.Array:
    tmp = data
    for field_tmp in field.split("/"):
        tmp = tmp[field_tmp]
    return tmp


def get_drift_time(
    data: ak.Array, end_time_field: str, t0_field: str | None
) -> ak.Array:
    """Get the drift time."""
    end_time = _get_nested_field(data, end_time_field)

    if t0_field is not None:
        drift_time = end_time - _get_nested_field(data, t0_field)
    else:
        drift_time = end_time

    return drift_time


def _select_data_in_slice(
    drift_time: ak.Array,
    energy: ak.Array,
    drift_slice: Slice,
) -> ak.Array:
    """Filter single-detector event data to one energy-drift-time slice."""
    return (
        ak.all(energy >= drift_slice.energy_range[0], axis=-1)
        & ak.all(energy <= drift_slice.energy_range[1], axis=-1)
        & ak.all(drift_time >= drift_slice.drift_time_range[0], axis=-1)
        & ak.all(drift_time <= drift_slice.drift_time_range[1], axis=-1)
    )


def lookup_wfs_indices(
    slices: list[Slice],
    *,
    evt_files: list[str],
    n_target: int,
    detector: str,
    t0_field: str | None = "spms/event_t0",
    end_time_field: str = "geds/psd/low_aoe/time",
) -> list[AttrsDict]:
    """Extract the indices of the waveforms to use in superpulse construction.

    The accumulation is stopped when either every slice has
    more than `n_target` waveforms of the average is more than
    3 times `n_target`.

    Parameters
    ----------
    slices
        A list of slices to extract wf indices for.
    evt_files
        List of evt files.
    n_target
        The maximum number of waveforms to select.
    detector
        The detector to use.
    t0_field
        Field for the start time, if `None` will be set to 0.
    end_time_field
        Field for the end-time of the drift time calculation.

    Returns
    -------
    a tuple of:
        a list with the same length as "slices", each element is an `AttrsDict` with three fields
        - "file_idx": the indice of the file lists containing the selected waveforms,
        - "hit_idx": the row of the files,
        - "n_sel": the number of selected waveforms.
        and a list of all drift times (of considered files)
    """
    output = [AttrsDict({"file_idx": [], "hit_idx": [], "n_sel": 0}) for _ in slices]

    all_drift_times = None
    for file_idx, evt_file in enumerate(evt_files):
        # early break to speed up
        m = np.mean([out.n_sel for out in output])
        if all(out.n_sel >= n_target for out in output) or (m >= 3 * n_target):
            break

        if file_idx % 100 == 0:
            msg = f"Reading file {file_idx} out of {len(evt_files)} {m} target ({n_target})"
            log.info(msg)

        evts = _read_and_sel_evts(evt_file, detector=detector, t0_field=t0_field)
        drift_time = get_drift_time(evts, end_time_field, t0_field)

        all_drift_times = (
            np.concatenate((all_drift_times, ak.flatten(drift_time)))
            if all_drift_times is not None
            else ak.flatten(drift_time).to_numpy()
        )

        for out_tmp, drift_slice in zip(output, slices, strict=True):
            if out_tmp.n_sel >= n_target:
                continue

            # evts in our slice
            evts_slice = evts[
                _select_data_in_slice(
                    drift_time, evts.geds.energy, drift_slice=drift_slice
                )
            ]
            hit_indices = ak.flatten(evts_slice.geds.hit_idx).to_list()

            out_tmp.hit_idx.extend(hit_indices)
            out_tmp.file_idx.extend(list(np.full_like(hit_indices, file_idx)))
            out_tmp.n_sel += len(hit_indices)

    return [
        AttrsDict(
            {
                "file_idx": out.file_idx[:n_target],
                "hit_idx": out.hit_idx[:n_target],
                "n_sel": len(out.hit_idx[:n_target]),
            }
        )
        for out in output
    ], (all_drift_times if all_drift_times is not None else np.array([], dtype=float))


def _get_dsp_config(dsp_config: str | dict | Path) -> dict:
    if isinstance(dsp_config, dict):
        config_dict = dsp_config
    elif isinstance(dsp_config, (str, Path)):
        config_dict = dbetto.utils.load_dict(dsp_config)
    else:
        msg = "dsp_config must be a dict or a path to a file"
        raise ValueError(msg)

    # modify the current so that the output is not cut
    config_dict["processors"]["curr"] = {
        "description": "differentiate waveform",
        "function": "avg_current",
        "module": "dspeed.processors",
        "args": [
            "wf_pz_win",
            "1",
            "curr(shape=len(wf_pz_win)-1,  period=wf_pz_win.period,offset=wf_pz_win.offset)",
        ],
        "unit": "ADC/sample",
    }
    return config_dict


def get_wfs_for_slice(
    raw_files: list[str],
    lh5_group: str,
    hit_indices: list[int],
    file_indices: list[int],
    *,
    dsp_config: str | dict | Path,
    charge_output: str = "wf_pz_win",
    current_output: str = "curr_av",
    bl_output: str = "bl_std_win",
    energy_output: str = "cuspEmax",
    align: str = "tp_aoe_max",
) -> AttrsDict | None:
    """Extract aligned charge and current waveforms for a set of slice events.

    Uses the DSP processing chain via ``WaveformBrowser``. After alignment each
    event has the same sampling rate and number of samples but a different
    x-offset. This function collects all events, finds the overlapping time
    region, and trims every waveform to that common window.

    Parameters
    ----------
    raw_files
        List of raw files.
    hit_indices
        The list of rows in each file to read.
    file_indices
        List of file indices to read.
    lh5_group
        HDF5 group containing the waveform table, e.g. ``"ch1084803/raw"``.
    dsp_config
        Path to the production DSP configuration JSON file.
    charge_output
        DSP output name for the charge waveform.
    current_output
        DSP output name for the current waveform.
    bl_output
        DSP output name for the baseline standard deviation.
    energy_output
        DSP output name for the energy,
    align
        DSP parameter used to align waveforms on the time axis.

    Returns
    -------
    A dictionary with the following fields (or None if no valid waveforms found):
        charge_times
            Common time axis for charge waveforms, shape ``(n_common_charge,)``.
        current_times
            Common time axis for current waveforms, shape ``(n_common_current,)``.
        charge_wfs
            Trimmed charge waveforms, shape ``(n_valid, n_common_charge)``.
        current_wfs
            Trimmed current waveforms, shape ``(n_valid, n_common_current)``.
        bl_std
            Baseline standard deviation for each event, shape ``(n_valid,)``.
        energy
            Energy estimator for each event, shape ``(n_valid,)``.
    """
    from dspeed.vis import WaveformBrowser  # noqa: PLC0415

    waveforms = []
    for file_idx in np.unique(file_indices):
        indices = [
            int(idx)
            for idx in np.array(hit_indices)[np.array(file_indices) == file_idx]
        ]

        if len(indices) == 0:
            continue

        dsp = _get_dsp_config(dsp_config)

        browser = WaveformBrowser(
            str(raw_files[file_idx]),
            lh5_group,
            dsp_config=dsp,
            lines=[charge_output, current_output, bl_output, energy_output],
            align=align,
        )

        browser.find_entry(indices, append=False)

        charge_lines = browser.lines.get(charge_output, [])
        current_lines = browser.lines.get(current_output, [])

        bl_std_vals = browser.lines.get(bl_output, [])
        energy_vals = browser.lines.get(energy_output, [])

        # First pass: collect valid events with their x- and y-data
        for i, (cl, il, el) in enumerate(
            zip(charge_lines, current_lines, energy_vals, strict=True)
        ):
            charge_y = cl.get_ydata() / float(el.get_ydata()[0])
            current_y = il.get_ydata() / float(el.get_ydata()[0])

            if np.any(np.isnan(charge_y)) or np.any(np.isnan(current_y)):
                log.debug("event %d: NaN in waveform, skipping", i)
                continue

            # alignment failed (skipped)
            if align is not None and (
                cl.get_xdata()[0] == 0.0 or il.get_xdata()[0] == 0.0
            ):
                continue

            entry = {
                "charge_t": cl.get_xdata(),
                "charge_y": charge_y,
                "current_t": il.get_xdata(),
                "current_y": current_y,
            }

            entry["bl_std"] = float(bl_std_vals[i].get_ydata()[0])
            entry["energy"] = float(el.get_ydata()[0])

            waveforms.append(entry)

    if len(waveforms) == 0:
        return None

    # Find the common time window (max of left edges, min of right edges)
    charge_t_min = max(-2000, *(e["charge_t"][0] for e in waveforms))
    charge_t_max = min(3000.0, *(e["charge_t"][-1] for e in waveforms))

    current_t_min = max(-2000, *(e["current_t"][0] for e in waveforms))
    current_t_max = min(3000.0, *(e["current_t"][-1] for e in waveforms))

    log.info("Current range across events: [%f, %f] ns", current_t_min, current_t_max)
    log.info("Charge range across events: [%f, %f] ns", charge_t_min, charge_t_max)

    if (current_t_max <= current_t_min) or (charge_t_max <= charge_t_min):
        log.warning(
            "No overlapping time window found for waveforms (current: [%f, %f], charge: [%f, %f]), skipping",
            current_t_min,
            current_t_max,
            charge_t_min,
            charge_t_max,
        )
        return None

    charge_times = np.arange(charge_t_min, charge_t_max)

    # Trim each event to the common window
    charge_list = []
    curr_list = []
    bl_std_list = []
    energy_list = []

    for e in waveforms:
        # interpolate charge to 1 ns
        charge_wf = np.interp(charge_times, e["charge_t"], e["charge_y"])

        i_mask = (e["current_t"] > current_t_min) & (e["current_t"] < current_t_max)

        charge_list.append(charge_wf)
        curr_list.append(np.array(e["current_y"][i_mask]))

        bl_std_list.append(e["bl_std"])
        energy_list.append(e["energy"])

    # Common time axes from the first event's trimmed region
    e0 = waveforms[0]

    current_times = e0["current_t"][
        (e0["current_t"] > current_t_min) & (e0["current_t"] < current_t_max)
    ]

    charge_stack = np.vstack(charge_list)
    curr_stack = np.vstack(curr_list)

    if len(charge_times) != np.shape(charge_stack)[1]:
        msg = (
            f"charge time axis length ({len(charge_times)}) does not match "
            f"number of samples in charge waveforms ({np.shape(charge_stack)[1]})."
        )
        raise ValueError(msg)

    if len(current_times) != np.shape(curr_stack)[1]:
        msg = (
            f"current time axis length ({len(current_times)}) does not match "
            f"number of samples in current waveforms ({np.shape(curr_stack)[1]})."
        )
        raise ValueError(msg)

    energy = np.array(energy_list)
    bl_std = np.array(bl_std_list)

    return AttrsDict(
        {
            "charge_times": charge_times,
            "current_times": current_times,
            "charge_wfs": charge_stack,
            "current_wfs": curr_stack,
            "bl_std": bl_std,
            "energy": energy,
        }
    )


def compute_chi2(
    charge_wfs: np.ndarray,
    superpulse: Superpulse,
    bl_std: np.ndarray,
    cuspEmax: np.ndarray,
) -> np.ndarray:
    """
    Compute the reduced chi-squared of each charge waveform vs the superpulse.

    The normalised noise is ``sigma_wf = bl_std / cuspEmax`` for each event, and the superpulse
      sigma is the mean of these values.

    Parameters
    ----------
    charge_wfs
        2D array of shape ``(n_events, n_samples)``.
    superpulse
        Preliminary superpulse.
    bl_std
        Per-event baseline standard deviation in ADC units, shape
        ``(n_events,)``.
    cuspEmax
        Per-event energy estimator in ADC units, shape ``(n_events,)``.

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

    # Normalised noise from DSP chain
    sigma_wfs = bl_std / cuspEmax  # (n_events,)
    sigma_sp = np.mean(sigma_wfs)

    # Vectorised chi2: residuals shape (n_events, n_samples), sigma_wfs broadcast from (n_events, 1)
    residuals_sq = (charge_wfs - sp_wf[np.newaxis, :]) ** 2
    variance = sigma_wfs[:, np.newaxis] ** 2 + sigma_sp**2
    chi2 = np.nansum(residuals_sq / variance, axis=1)

    return chi2 / n_samples


def write_superpulses(
    superpulses: dict[Slice, Superpulse],
    output_path: str,
    detector: str,
    *,
    wo_mode: str = "write_safe",
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
                drift_time_center     [Scalar, ns]
                drift_time_lo         [Scalar, ns]
                drift_time_hi         [Scalar, ns]
                energy_lo             [Scalar, keV]
                energy_hi             [Scalar, keV]
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
    wo_mode
        Write mode for lh5 files. Default: "write_safe".
    """
    for idx, (sl, sp) in enumerate(superpulses.items()):
        # Build group name: {detector}/dt_{lo}_{hi}_ns

        for v in sl.drift_time_range:
            value_float = float(v)
            if not np.isfinite(value_float) or not value_float.is_integer():
                msg = f"drift time bounds must be finite integers, got {sl.drift_time_range}"
                raise ValueError(msg)

        dt_lo = int(sl.drift_time_range[0])
        dt_hi = int(sl.drift_time_range[1])
        group = f"{detector}/dt_{dt_lo}_{dt_hi}_ns"

        lgdo_struct = sp.to_lgdo()
        lh5.write(
            lgdo_struct, group, output_path, wo_mode=wo_mode if idx == 0 else "append"
        )

        log.debug("wrote %s to %s", group, output_path)

    log.info(
        "wrote %d slices for %s to %s",
        len(superpulses),
        detector,
        output_path,
    )


def read_superpulses(
    path: str,
    detector: str,
) -> dict[Slice, Superpulse]:
    """Read superpulses written by :func:`write_superpulses`.

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
            energy_range=(
                float(group["energy_lo"].value),
                float(group["energy_hi"].value),
            ),
            drift_time_range=(
                float(group["drift_time_lo"].value),
                float(group["drift_time_hi"].value),
            ),
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
    *,
    xlims=(-1000, 3000),
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
        Selected waveform arrays after chi2 cut, shape ``(n_events, n_samples)``.
    superpulse
        The final superpulse to overlay.

    Returns
    -------
    fig
    (ax_charge, ax_current) axes for the two panels.
    """
    sl = superpulse.slice

    # Shared x range: intersection of the two time axes

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
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
    ax1.set_xlim(*xlims)
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
    ax2.set_xlim(*xlims)
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
    elif curve == "current":
        sp_times = final_superpulse.current_time_axis
        sp_wf = final_superpulse.current_wf
        ylabel = "d(ADC/cuspEmax)/dt"
    else:
        msg = f"invalid curve '{curve}'; expected one of ('charge', 'current')"
        raise ValueError(msg)

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
) -> tuple:
    """Plot superpulses from all drift-time slices, color-coded by drift time.

    Parameters
    ----------
    lh5_file
        Path to the LH5 file produced by :func:`write_superpulses`.
    detector
        Detector name (top-level group in the file).
    curve
        ``"charge"`` or ``"current"``.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax  : matplotlib.axes.Axes
    """
    import matplotlib.colors as mcolors  # noqa: PLC0415

    if curve not in ("charge", "current"):
        msg = f"curve must be 'charge' or 'current', got {curve!r}"
        raise ValueError(msg)

    groups = sorted(
        lh5.ls(lh5_file, f"{detector}/"),
        key=lambda g: int(g.split("_")[1]),
    )

    # collect all slices first so we can normalise the colormap to actual data
    slices = []
    for group in groups:
        struct = lh5.read(group, lh5_file)
        slices.append(
            {
                "times": struct[f"{curve}_time_axis"].nda,
                "wf": struct[f"{curve}_wf"].nda,
                "dt_center": struct["drift_time_center"].value,
                "dt_lo": struct["drift_time_lo"].value,
                "dt_hi": struct["drift_time_hi"].value,
                "e_lo": struct["energy_lo"].value,
                "e_hi": struct["energy_hi"].value,
                "n_events": struct["n_events_final"].value,
            }
        )

    dt_centers = [s["dt_center"] for s in slices]
    dt_min, dt_max = min(dt_centers), max(dt_centers)
    norm = mcolors.Normalize(vmin=dt_min, vmax=dt_min + (dt_max - dt_min) / 0.85)
    cmap = colormaps["viridis"]

    ylabel = "ADC / cuspEmax" if curve == "charge" else "d(ADC/cuspEmax)/dt"
    e_lo, e_hi = slices[0]["e_lo"], slices[0]["e_hi"]

    fig, ax = plt.subplots(figsize=(8, 6))

    for s in slices:
        ax.plot(s["times"], s["wf"], color=cmap(norm(s["dt_center"])), linewidth=2)
    ax.set_xlim(-2000, 1000)

    # 10% amplitude reference line (charge only)
    if curve == "charge":
        ax.axhline(0.10, color="gray", linestyle="--", linewidth=1)
        xlo, xhi = ax.get_xlim()
        x_ann = xlo + 0.95 * (xhi - xlo)
        ax.annotate(
            "10% amplitude",
            xy=(x_ann, 0.10),
            xytext=(x_ann, 0.107),
            fontsize=10,
            color="gray",
            va="bottom",
            ha="right",
        )
    ax.set_title(f"{detector} - Average {curve} pulses  |  E = [{e_lo}, {e_hi}] keV")
    ax.set_xlabel("Time [ns]")
    ax.set_ylabel(ylabel)
    ax.grid(True, color="grey", linestyle="--", alpha=0.2)

    # colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Drift Time [ns]", labelpad=10)

    # inset
    ax_inset = ax.inset_axes([0.08, 0.5, 0.45, 0.45])
    for s in slices:
        ax_inset.plot(
            s["times"],
            s["wf"],
            color=cmap(norm(s["dt_center"])),
            linewidth=1.8,
        )
    ax_inset.tick_params(labelsize=7)
    ax_inset.grid(True, color="grey", linestyle="--", alpha=0.2)

    if curve == "charge":
        ax_inset.set_xlim(-2000, -100)
        ax_inset.set_ylim(-0.001, 0.105)
        # 1% amplitude reference line
        ax_inset.axhline(0.01, color="gray", linestyle="--", linewidth=0.8)
        xlo, xhi = ax_inset.get_xlim()
        x_ann = xlo + 0.05 * abs(xhi - xlo)
        ax_inset.annotate(
            "1% amplitude", xy=(x_ann, 0.012), fontsize=8, color="gray", ha="left"
        )
    else:
        ax_inset.set_xlim(-50, 50)
        xlo, xhi = ax_inset.get_xlim()
        y_in_view = [
            v
            for s in slices
            for t, v in zip(s["times"], s["wf"], strict=True)
            if xlo <= t <= xhi
        ]
        if y_in_view:
            ax_inset.set_ylim(min(y_in_view) - 0.002, max(y_in_view) + 0.002)

    fig.tight_layout()
    return fig, ax


def plot_current_superpulses_fwhm_and_amplitude(
    lh5_file: str,
    detector: str,
    dt_range_tuning: tuple[float, float] | None = None,
) -> tuple:
    """Plot FWHM and peak amplitude of current superpulses vs drift time.

    Two-row figure sharing the x-axis (drift time center). Color coding
    matches :func:`plot_superpulses`.

    NOTE: dt_range_tuning shows the contiguous span of slices used
    # for the fit; individual points within this range may not all
    # have been used if max_num_superpulses further truncated them

    Parameters
    ----------
    lh5_file
        Path to the LH5 file produced by :func:`write_superpulses`.
    detector
        Detector name (top-level group in the file).
    dt_range_tuning
        Optional tuple of (dt_lo, dt_hi) to indicate the drift time range used for tuning.

    Returns
    -------
    fig : matplotlib.figure.Figure
    (ax_fwhm, ax_amp) : tuple of Axes
    """
    import matplotlib.colors as mcolors  # noqa: PLC0415

    groups = sorted(
        lh5.ls(lh5_file, f"{detector}/"),
        key=lambda g: int(g.split("_")[1]),
    )

    slices = []
    for group in groups:
        struct = lh5.read(group, lh5_file)
        times = struct["current_time_axis"].nda
        wf = struct["current_wf"].nda
        dt_center = struct["drift_time_center"].value

        # Current amplitude
        amax = np.nanmax(wf)
        if np.isnan(amax):
            log.warning(
                "slice %s: current waveform contains only NaNs, skipping", group
            )
            continue

        # FWHM via linear interpolation at half-max crossings
        half_max = amax / 2.0
        idx_peak = np.argmax(wf)
        shifted = wf - half_max

        sign_changes = np.where(np.diff(np.sign(shifted)))[0]

        left = sign_changes[sign_changes < idx_peak]
        right = sign_changes[sign_changes >= idx_peak]

        if len(left) == 0 or len(right) == 0:
            log.warning("slice %s: half-max crossing not found, skipping", group)
            continue

        il = left[-1]
        t_left = times[il] - shifted[il] * (times[il + 1] - times[il]) / (
            shifted[il + 1] - shifted[il]
        )
        ir = right[0]
        t_right = times[ir] - shifted[ir] * (times[ir + 1] - times[ir]) / (
            shifted[ir + 1] - shifted[ir]
        )

        fwhm = t_right - t_left

        slices.append(
            {
                "dt_center": dt_center,
                "amax": amax,
                "fwhm": fwhm,
                "e_lo": struct["energy_lo"].value,
                "e_hi": struct["energy_hi"].value,
            }
        )

    dt_centers = [s["dt_center"] for s in slices]
    dt_min, dt_max = min(dt_centers), max(dt_centers)
    norm = mcolors.Normalize(vmin=dt_min, vmax=dt_min + (dt_max - dt_min) / 0.85)
    cmap = colormaps["viridis"]

    e_lo, e_hi = slices[0]["e_lo"], slices[0]["e_hi"]

    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        figsize=(8, 4),
        sharex=True,
        gridspec_kw={"hspace": 0.08},
        constrained_layout=True,
    )

    for s in slices:
        color = cmap(norm(s["dt_center"]))
        ax1.plot(s["dt_center"], s["fwhm"], "o", color=color, markersize=7)
        ax2.plot(s["dt_center"], s["amax"], "o", color=color, markersize=7)

    ax1.set_ylabel("Current FWHM [ns]")
    ax1.grid(True, color="grey", linestyle="--", alpha=0.2)
    ax1.set_ylim(
        bottom=np.min([s["fwhm"] for s in slices]) - 2,
        top=np.max([s["fwhm"] for s in slices]) + 2,
    )

    ax2.set_ylabel("Current amplitude [A.U.]")
    ax2.set_xlabel("Drift time [ns]")
    ax2.grid(True, color="grey", linestyle="--", alpha=0.2)
    ax2.set_ylim(
        bottom=np.min([s["amax"] for s in slices]) - 0.002,
        top=np.max([s["amax"] for s in slices]) + 0.002,
    )

    fig.suptitle(
        f"{detector} - Response uniformity  | E = [{e_lo}, {e_hi}] keV",
        fontsize=13,
    )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=[ax1, ax2], pad=0.02)
    cbar.set_label("Drift Time [ns]", labelpad=10)

    if dt_range_tuning is not None:
        for ax in (ax1, ax2):
            ax.axvline(dt_range_tuning[0], color="black", linestyle="--", linewidth=0.8)
            ax.axvline(dt_range_tuning[1], color="black", linestyle="--", linewidth=0.8)
            xlim = ax.get_xlim()
            ax.axvspan(xlim[0], dt_range_tuning[0], color="grey", alpha=0.2)
            ax.axvspan(dt_range_tuning[1], xlim[1], color="grey", alpha=0.2)
            ax.set_xlim(xlim)

    return fig, (ax1, ax2)
