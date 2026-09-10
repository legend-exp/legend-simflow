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
from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from pathlib import Path

import awkward as ak
import lh5
import numpy as np
import pyg4ometry
import pygeomtools
import reboost.hpge
import reboost.math
import reboost.units
from numpy.typing import ArrayLike

from legendsimflow import nersc, utils

from . import patterns
from .utils import SimflowConfig

log = logging.getLogger(__name__)


def get_senstables(
    geom: pyg4ometry.geant4.Registry, det_type: str | None = None
) -> list[str]:
    sensvols = pygeomtools.detectors.get_all_senstables(geom)
    if det_type is not None:
        return [k for k, v in sensvols.items() if v.detector_type == det_type]
    return list(sensvols.keys())


def load_hpge_realistic_psl(
    config: SimflowConfig,
    det_name: str,
    runid: str,
    waveform_angles: Sequence[int] = (0,),
) -> (
    tuple[
        dict[int, reboost.hpge.HPGeRZField],
        dict[int, reboost.hpge.HPGePulseShapeLibrary],
    ]
    | None
):
    """Load HPGe PSL from disk.

    Loads the waveforms as well as drift-time maps for
    both coordinates.

    Parameters
    ----------
    config
        Simflow configuration object.
    det_name
        HPGe detector name.
    runid
        Run identifier.
    waveform_angles
        Crystal-axis angles, in degrees, for which to load the (large) waveform
        pulse-shape libraries. Defaults to the 0-degree axis only, which is
        the single axis consumed by :func:`extract_detailed_psd_observables`.
        The drift-time maps are always loaded for both axes regardless.

    Returns
    -------
    A tuple of two dictionaries: the first contains the drift-time maps for different crystal axes, keyed by angle.
    The second contains the corresponding pulse-shape libraries.
    If no valid maps are found, ``None`` is returned for both.

    """
    psl_file = nersc.dvs_ro(
        config,
        patterns.output_realistic_psl_merged_filename(
            config,
            runid=runid,
        ),
    )

    if (
        Path(psl_file).exists()
        and len(lh5.ls(psl_file, f"{det_name}/drift_time_*")) >= 2
    ):
        log.debug("loading drift-time maps from %s", psl_file)
        # both axes needed: drift_time_crystal_axes() blends them (and small)
        dt_map = reboost.hpge.load_hpge_drift_time_maps(
            psl_file, det_name, bounds_error=False
        )
        # waveforms dominate memory (PSL file is O(10 GB)); load only the axes
        # used downstream (extract_detailed_psd_observables() is single-axis).
        # float32 is enough for ~1% A/E; halves the templates if read as f64
        waveform_map = reboost.hpge.load_hpge_pulse_shape_libraries(
            psl_file,
            det_name,
            angles=waveform_angles,
            out_of_bounds_val=0,
            dtype=np.float32,
        )

    else:
        msg = (
            f"no valid time maps found for {det_name} in {psl_file}, "
            "drift time will be set to NaN"
        )
        log.warning(msg)
        dt_map = None
        waveform_map = None

    return dt_map, waveform_map


def load_hpge_dtmaps(
    config: SimflowConfig, det_name: str, runid: str
) -> dict[int, reboost.hpge.HPGeRZField] | None:
    """Load HPGe drift-time maps from disk.

    Automatically finds and loads drift-time maps for crystal axes <100> <110>.
    If no map is found, ``None`` is returned.

    Parameters
    ----------
    config
        Simflow configuration object.
    det_name
        HPGe detector name.
    runid
        Run identifier.

    """
    hpge_dtmap_file = patterns.output_dtmap_merged_filename(
        config,
        runid=runid,
    )

    if len(lh5.ls(hpge_dtmap_file, f"{det_name}/drift_time_*")) >= 2:
        log.debug("loading drift-time maps")
        dt_map = reboost.hpge.load_hpge_drift_time_maps(
            hpge_dtmap_file, det_name, bounds_error=False
        )
    else:
        msg = (
            f"no valid time maps found for {det_name} in {hpge_dtmap_file}, "
            "drift time will be set to NaN"
        )
        log.warning(msg)
        dt_map = None

    return dt_map


def extract_psd_observables(
    chunk: ak.Array,
    edep_active: ak.Array,
    energy: ak.Array,
    dt_map: dict[int, reboost.hpge.HPGeRZField],
    currmod_pars: Mapping,
    det_loc: pyg4ometry.gdml.Defines.Position,
    *,
    aoe_res: ArrayLike,
    aoe_mean: ArrayLike,
    psdcuts: Mapping,
    current_reso: float,
) -> ak.Array:
    """Extract PSD observables for a chunk of events in an HPGe detector.

    This function calculates the A/E observable, its classifier, and the single-site flag for a chunk of events in an HPGe detector, using the provided drift-time maps and current model parameters.

    Parameters
    ----------
    chunk
        Awkward array containing the events to process. Must have fields 'xloc', 'yloc', 'zloc'.
    edep_active
        Energy deposited in the active volume per hit, used for A/E calculation.
    energy
        Energy deposited in the active volume, used for A/E calculation.
    dt_map
        Dictionary of drift-time maps for different crystal axes, as returned by `load_hpge_dtmaps()`.
    currmod_pars
        Dictionary of parameters for the current model, (see
        :func:`reboost.hpge.get_current_template`)
    det_loc
        Position of the detector in the global coordinate system, used for drift time correction.
    det_name
        Name of the detector.
    aoe_res
        A/E resolution (sigma) used for A/E classifier calculation, typically determined from data.
    aoe_mean
        A/E mean used in the A/E classifier calculation, typically from fitting simulated data.
    psdcuts
        Dictionary containing the low and high side cuts for the A/E classifier to determine single-site events.
    current_reso
        Standard deviation of the Gaussian noise to smear the maximum current, representing the current resolution of the detector.


    Returns
    -------
    an `ak.Array` with fields:
        - aoe: A/E observable for each event.
        - aoe_class: A/E classifier (normalized to resolution) for each event.
        - is_single_site: boolean flag indicating whether the event is classified as single-site based on A/E cuts.
        - t_max: drift time at the position of maximum current, useful for further PSD or analysis.
    """
    _drift_time = reboost.hpge.drift_time_crystal_axes(
        chunk.xloc, chunk.yloc, chunk.zloc, dt_map, coord_offset=det_loc
    )
    utils.check_nans_leq(_drift_time, "_drift_time", 0.01, min_entries=1000)

    _a_max_true = hpge_max_current(edep_active, _drift_time, currmod_pars)

    utils.check_nans_leq(_a_max_true, "_a_max_true", 0.01, min_entries=1000)

    # Apply current resolution smearing based on configured A/E noise parameters
    _a_max = gauss_smear(_a_max_true, current_reso)

    # finally calculate A/E, comparable to the A/E in data
    # corrected for energy dependence
    aoe = _a_max / energy

    aoe_corr = aoe - aoe_mean + 1

    # ...and A/E classifier
    # NOTE: we use the resolution determined from data here instead
    # of the intrinsic simulated ones due to noise
    aoe_class = (aoe_corr - 1) / aoe_res

    # ...and PSD flag
    is_single_site = aoe_class > psdcuts.aoe.low_side

    # also calculate drift time at A position
    # FIXME: this is wasting compute resources, max_current should
    # return (maxA, t_maxA)
    t_max = hpge_max_current(
        edep_active,
        _drift_time,
        currmod_pars,
        return_mode="max_time",
    )

    return ak.Array(
        {
            "aoe": aoe,
            "aoe_corr": aoe_corr,
            "aoe_class": aoe_class,
            "is_single_site": is_single_site,
            "t_max": t_max,
        }
    )


def extract_detailed_psd_observables(
    chunk: ak.Array,
    edep_active: ak.Array,
    energy: ak.Array,
    dt_map: dict[int, reboost.hpge.HPGeRZField],
    pulse_shape_lib: reboost.hpge.HPGePulseShapeLibrary,
    det_loc: pyg4ometry.gdml.Defines.Position,
    *,
    aoe_res: ArrayLike,
    aoe_mean: ArrayLike,
    psdcuts: Mapping,
    current_reso: float | None = None,
) -> ak.Array:
    """Extract PSD observables for a chunk of events in an HPGe detector.

    This function calculates the A/E observable, its classifier, and the single-site flag for a chunk of events in an HPGe detector, using the provided drift-time maps and current model parameters.

    Parameters
    ----------
    chunk
        Awkward array containing the events to process. Must have fields 'xloc', 'yloc', 'zloc'.
    edep_active
        Energy deposited in the active volume per hit, used for A/E calculation.
    energy
        Energy deposited in the active volume, used for A/E calculation.
    dt_map
        Dictionary of drift-time maps for different crystal axes, as returned by `load_hpge_dtmaps()`.
    pulse_shape_lib
        Dictionary of waveform templates for different crystal axes, as returned by `load_hpge_realistic_psl()`.
    det_loc
        Position of the detector in the global coordinate system, used for drift time correction.
    det_name
        Name of the detector.
    aoe_res
        A/E resolution (sigma) used for A/E classifier calculation, typically determined from data.
    aoe_mean
        A/E mean used in the A/E classifier calculation, typically from fitting simulated data.
    psdcuts
        Dictionary containing the low and high side cuts for the A/E classifier to determine single-site events.
    current_reso
        Standard deviation of the Gaussian noise to smear the maximum current, representing the current resolution of the detector.

    Returns
    -------
    an `ak.Array` with fields:
        - aoe: A/E observable for each event.
        - aoe_class: A/E classifier (normalized to resolution) for each event.
        - is_single_site: boolean flag indicating whether the event is classified as single-site based on A/E cuts.
        - t_max: drift time at the position of maximum current, useful for further PSD or analysis.
    """
    # Convert det_loc to pint Quantity
    det_loc_pint = reboost.units.pg4_to_pint(det_loc)

    # Use reboost.units to get conversion factors for chunk coordinates
    # This handles the case when chunk has units attached (with_units=True)
    xloc_conv = reboost.units.units_convfact(chunk.xloc, det_loc_pint.units)
    yloc_conv = reboost.units.units_convfact(chunk.yloc, det_loc_pint.units)
    zloc_conv = reboost.units.units_convfact(chunk.zloc, det_loc_pint.units)

    # Unwrap LGDO/pint if present
    xloc, _ = reboost.units.unwrap_lgdo(chunk.xloc)
    yloc, _ = reboost.units.unwrap_lgdo(chunk.yloc)
    zloc, _ = reboost.units.unwrap_lgdo(chunk.zloc)

    _x = xloc * xloc_conv - det_loc_pint[0].m
    _y = yloc * yloc_conv - det_loc_pint[1].m

    _z = reboost.units.attach_units(1000 * (zloc * zloc_conv - det_loc_pint[2].m), "mm")
    _r = reboost.units.attach_units(1000 * np.sqrt(_x**2 + _y**2), "mm")

    _drift_time = reboost.hpge.drift_time_crystal_axes(
        chunk.xloc, chunk.yloc, chunk.zloc, dt_map, coord_offset=det_loc
    )

    utils.check_nans_leq(_drift_time, "_drift_time", 0.1, min_entries=1000)
    _a_max_true = reboost.hpge.maximum_current(
        edep_active,
        _drift_time,
        times=None,
        r=_r,
        z=_z,
        template=pulse_shape_lib,
        return_mode="current",
    )

    utils.check_nans_leq(_a_max_true, "_a_max_true", 0.01, min_entries=1000)

    # Apply current resolution smearing based on configured A/E noise parameters
    _a_max = gauss_smear(_a_max_true, current_reso)

    # finally calculate A/E, comparable to the A/E in data
    # corrected for energy dependence
    aoe = _a_max / energy

    aoe_corr = aoe - aoe_mean + 1
    # ...and A/E classifier
    # NOTE: we use the resolution determined from data here instead
    # of the intrinsic simulated ones due to noise
    aoe_class = (aoe_corr - 1) / aoe_res

    # ...and PSD flag
    is_single_site = aoe_class > psdcuts.aoe.low_side
    is_high_aoe = aoe_class > psdcuts.aoe.high_side
    is_bb_like = is_single_site & (~is_high_aoe)

    # also calculate drift time at A position
    # FIXME: this is wasting compute resources, max_current should
    # return (maxA, t_maxA)
    t_max = reboost.hpge.maximum_current(
        edep_active,
        _drift_time,
        times=None,
        r=_r,
        z=_z,
        template=pulse_shape_lib,
        return_mode="max_time",
    )

    return ak.Array(
        {
            "aoe": aoe,
            "aoe_corr": aoe_corr,
            "aoe_class": aoe_class,
            "is_single_site": is_single_site,
            "is_bb_like": is_bb_like,
            "is_high_aoe": is_high_aoe,
            "t_max": t_max,
        }
    )


def hpge_max_current(
    edep: ak.Array,
    drift_time: ak.Array,
    currmod_pars: Mapping,
    **kwargs,
) -> ak.Array:
    """Calculate the maximum of the current pulse.

    Parameters
    ----------
    edep
        energy deposited at each step.
    drift_time
        drift time of each energy deposit.
    currmod_pars
        dictionary storing the parameters of the current model (see
        :func:`reboost.hpge.get_current_template`)
    kwargs
        forwarded to :func:`reboost.hpge.maximum_current`.

    """
    # current pulse template domain in ns (step is 1 ns)
    t_domain = {"low": -1000, "high": 4000, "step": 1}

    # instantiate the template
    a_tmpl, times = reboost.hpge.get_current_template(
        **t_domain,
        mean_aoe=1,  # set the maximum of the template to unity, so the A/E will be calibrated
        **currmod_pars,
    )
    # and calculate the maximum current
    return reboost.hpge.maximum_current(
        edep,
        drift_time,
        template=a_tmpl,
        times=times,
        **kwargs,
    )


def gauss_smear(arr_true: ak.Array, arr_reso: ak.Array) -> ak.Array:
    """Smear values with expected resolution.

    Samples from gaussian and shifts negative values to a fixed, tiny positive
    value.
    """
    arr_smear = reboost.math.gaussian_sample(
        arr_true,
        arr_reso,
    )

    # energy can't be negative as a result of smearing
    return ak.where((arr_smear <= 0) & (arr_true >= 0), np.finfo(float).tiny, arr_smear)
