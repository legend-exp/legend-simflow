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
from typing import NamedTuple

import awkward as ak
import lh5
import numpy as np
import pint
import pyg4ometry
import pygeomhpges
import pygeomtools
import reboost.hpge
import reboost.math
import reboost.units
from dbetto import AttrsDict
from dbetto.utils import load_dict
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


def read_detector_origins(stp_file: str | Path) -> dict[str, pint.Quantity]:
    """Read the detector origins stored in a remage output file.

    Returns a mapping of detector name to its origin ``[x, y, z]`` in the
    global coordinate system, as a :class:`pint.Quantity` in metres.

    Parameters
    ----------
    stp_file
        Path to a remage output file.
    """
    # FIXME: units should be already present, to be fixed in remage
    u = pint.UnitRegistry()
    det_loc = lh5.read("detector_origins", stp_file)
    return {
        k: [v[field].value for field in ("xloc", "yloc", "zloc")] * u.m
        for k, v in det_loc.items()
    }


def hpge_active_energy(
    chunk: ak.Array,
    pyobj: pygeomhpges.HPGe,
    det_loc: pint.Quantity,
    fccd_in_mm: float,
    dead_layer_fraction: float,
) -> tuple[ak.Array, ak.Array]:
    """Apply the HPGe active-volume model to a chunk of remage steps.

    Computes the distance of each step to the n+ surface of the detector (see
    :func:`reboost.hpge.distance_to_surface`), the corresponding
    charge-collection efficiency with a piecewise-linear dead-layer model (see
    :func:`reboost.math.piecewise_linear_activeness`) and weights
    the deposited energies accordingly.

    Parameters
    ----------
    chunk
        Awkward array of remage steps. Must have fields ``xloc``, ``yloc``,
        ``zloc``, ``edep`` and ``dist_to_surf``.
    pyobj
        The HPGe geometry object.
    det_loc
        Origin of the detector in the global coordinate system.
    fccd_in_mm
        Full charge-collection depth of the detector, in mm.
    dead_layer_fraction
        Fraction of the dead layer thickness at which the linear ramp in
        charge-collection efficiency starts.

    Returns
    -------
    edep_active
        Per-step energy deposits weighted by the activeness.
    energy_true
        Per-event active energy, sum of `edep_active` over the steps.
    """
    distance_to_nplus = reboost.hpge.distance_to_surface(
        chunk.xloc,
        chunk.yloc,
        chunk.zloc,
        pyobj,
        det_loc,
        distances_precompute=chunk.dist_to_surf,
        precompute_cutoff=(fccd_in_mm + 1),
        surface_type="nplus",
    )

    activeness = reboost.math.piecewise_linear_activeness(
        distance_to_nplus,
        fccd_in_mm=fccd_in_mm,
        dlf=dead_layer_fraction,
    )

    edep_active = chunk.edep * activeness
    return edep_active, ak.sum(edep_active, axis=-1)


class HPGePSDInputs(NamedTuple):
    """The per-detector, per-run inputs of the HPGe PSD simulation.

    See :func:`load_hpge_psd_inputs`.
    """

    dt_map: dict[int, reboost.hpge.HPGeRZField] | None
    """Drift-time maps per crystal axis (single-template method)."""
    currmod_pars: Mapping | None
    """Current-pulse model parameters (single-template method)."""
    current_reso: float | None
    """Standard deviation of the Gaussian noise smearing the maximum current."""
    psl_dt_maps: dict[int, reboost.hpge.HPGeRZField] | None
    """Drift-time maps per crystal axis of the realistic pulse-shape library."""
    realistic_psl: dict[int, reboost.hpge.HPGePulseShapeLibrary] | None
    """Realistic pulse-shape library waveforms per crystal axis."""

    @property
    def can_model_psd(self) -> bool:
        """Whether the single-template A/E can be simulated."""
        return (
            self.dt_map is not None
            and self.currmod_pars is not None
            and self.current_reso is not None
        )

    @property
    def can_model_psd_with_psl(self) -> bool:
        """Whether the pulse-shape-library A/E can be simulated."""
        return self.realistic_psl is not None and self.current_reso is not None


def load_hpge_psd_inputs(
    config: SimflowConfig,
    det_name: str,
    runid: str,
    *,
    simulate_psd: bool,
    simulate_psd_with_psl: bool,
    currmod_pars_all: Mapping | None = None,
) -> HPGePSDInputs:
    """Load the inputs of the HPGe PSD simulation for a detector in a run.

    Reads (from the known `par` step file patterns) the drift-time maps, the
    current-pulse model parameters and the realistic pulse-shape library of
    `det_name` in `runid`, as requested by the hit-tier settings. Missing
    inputs are returned as ``None``; inspect
    :attr:`HPGePSDInputs.can_model_psd` and
    :attr:`HPGePSDInputs.can_model_psd_with_psl` to know which PSD method can
    be simulated.

    Parameters
    ----------
    config
        Simflow configuration object.
    det_name
        HPGe detector name.
    runid
        Run identifier.
    simulate_psd
        Whether the single-template PSD simulation is enabled.
    simulate_psd_with_psl
        Whether the pulse-shape-library PSD simulation is enabled.
    currmod_pars_all
        Pre-loaded content of the merged current-pulse model file of `runid`
        (see :func:`legendsimflow.patterns.output_currmod_merged_filename`).
        Loaded from disk if ``None``.
    """
    dt_map = load_hpge_dtmaps(config, det_name, runid) if simulate_psd else None

    psl_dt_maps = realistic_psl = None
    if simulate_psd_with_psl:
        psl_dt_maps, realistic_psl = load_hpge_realistic_psl(config, det_name, runid)

    # the current-pulse model (noise smearing) is a product of the
    # single-template PSD chain, needed by both methods
    currmod_pars = current_reso = None
    if simulate_psd:
        if currmod_pars_all is None:
            currmod_pars_all = AttrsDict(
                load_dict(
                    nersc.dvs_ro(
                        config,
                        patterns.output_currmod_merged_filename(config, runid=runid),
                    )
                )
            )
        pars = currmod_pars_all.get(det_name, None)
        if pars is not None:
            currmod_pars = pars.get("current_pulse_pars", None)
            current_reso = pars["current_reso"] / pars["mean_aoe"]
        else:
            log.warning("no current-pulse model found for %s in %s", det_name, runid)

    if simulate_psd_with_psl and realistic_psl is not None and current_reso is None:
        log.warning(
            "the pulse-shape-library A/E of %s in %s cannot be simulated without "
            "the current-pulse model (noise smearing), enable simulate_psd",
            det_name,
            runid,
        )

    return HPGePSDInputs(dt_map, currmod_pars, current_reso, psl_dt_maps, realistic_psl)


def compute_hpge_psd_observables(
    chunk: ak.Array,
    edep_active: ak.Array,
    energy: ak.Array,
    det_loc: pint.Quantity,
    inputs: HPGePSDInputs,
    *,
    aoe_res: ArrayLike,
    aoe_mean: ArrayLike,
    aoe_mean_psl: ArrayLike,
    psdcuts: Mapping | None = None,
) -> tuple[ak.Array, ak.Array]:
    """Compute the HPGe PSD observables of a chunk of events with all enabled methods.

    Dispatches to :func:`extract_psd_observables` (single-template method) and
    :func:`extract_detailed_psd_observables` (pulse-shape-library method)
    according to what `inputs` allow to simulate. The observables of a method
    that cannot be simulated are filled with NaN (or ``False`` for flags).

    Parameters
    ----------
    chunk
        Awkward array of the events to process.
    edep_active
        Per-step energy deposits weighted by the activeness.
    energy
        Per-event energy used to compute A/E.
    det_loc
        Origin of the detector in the global coordinate system.
    inputs
        The PSD inputs of the detector in the run, see
        :func:`load_hpge_psd_inputs`.
    aoe_res
        A/E resolution (sigma) used for the A/E classifier.
    aoe_mean
        A/E mean energy-dependence model evaluated on the events
        (single-template method).
    aoe_mean_psl
        As `aoe_mean`, for the pulse-shape-library method.
    psdcuts
        Low and high side cuts of the A/E classifier. When omitted both are
        taken as zero, so the returned single-site flags are meaningless: only
        pass it if they are of interest.

    Returns
    -------
    psd_fields, psd_fields_detailed
        The observables of the single-template and pulse-shape-library methods
        (see the dispatched functions for the fields).
    """
    if psdcuts is None:
        psdcuts = AttrsDict({"aoe": {"low_side": 0.0, "high_side": 0.0}})

    n = len(chunk)
    psd_fields = ak.Array(
        {
            "aoe": np.full(n, np.nan),
            "aoe_class": np.full(n, np.nan),
            "aoe_corr": np.full(n, np.nan),
            "is_single_site": np.full(n, False),
            "t_max": np.full(n, np.nan),
        }
    )
    psd_fields_detailed = ak.Array(
        {
            "aoe": np.full(n, np.nan),
            "aoe_class": np.full(n, np.nan),
            "aoe_corr": np.full(n, np.nan),
            "is_single_site": np.full(n, False),
            "is_bb_like": np.full(n, False),
            "is_high_aoe": np.full(n, False),
            "t_max": np.full(n, np.nan),
        }
    )

    if inputs.can_model_psd:
        psd_fields = extract_psd_observables(
            chunk,
            edep_active,
            energy,
            inputs.dt_map,
            inputs.currmod_pars,
            det_loc,
            aoe_res=aoe_res,
            aoe_mean=aoe_mean,
            psdcuts=psdcuts,
            current_reso=inputs.current_reso,
        )

    if inputs.can_model_psd_with_psl:
        psd_fields_detailed = extract_detailed_psd_observables(
            chunk,
            edep_active,
            energy,
            inputs.psl_dt_maps,
            inputs.realistic_psl[0],
            det_loc,
            aoe_res=aoe_res,
            aoe_mean=aoe_mean_psl,
            psdcuts=psdcuts,
            current_reso=inputs.current_reso,
        )

    return psd_fields, psd_fields_detailed


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

    # the coordinates are now in the units of the detector origin, whatever
    # those are: label them so that reboost converts them itself
    det_loc_units = str(det_loc_pint.units)
    _z = reboost.units.attach_units(zloc * zloc_conv - det_loc_pint[2].m, det_loc_units)
    _r = reboost.units.attach_units(np.sqrt(_x**2 + _y**2), det_loc_units)

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
