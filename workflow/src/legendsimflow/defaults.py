# Copyright (C) 2026 Luigi Pertoldi <gipert@pm.me>
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
"""Default values of the optional Simflow configuration fields.

A field listed here is optional: a value set by the user replaces the default.
Nested blocks are merged field by field, lists are replaced as a whole. A field
not listed here is required, unless its documentation says otherwise.

The defaults are filled in when the configuration is read, by
:func:`.utils.init_simflow_context`, :func:`.metadata.get_tier_settings` and
:func:`.metadata.get_par_settings`. The meaning of each field is documented in
the user manual.
"""

from __future__ import annotations

import copy
from collections.abc import Mapping

from dbetto import AttrsDict

#: Simflow configuration file (``simflow-config.yaml``).
CONFIG: dict = {
    "simlist": "all",
    "simlist_exclude": [],
    "runlist": None,
    "make_steps": ["vtx", "stp", "par", "opt", "hit", "evt", "cvt", "pdf"],
    "metadata_version": None,
    "validate_l200data": True,
    "simflow_rng_seed": 0,
    "benchmark": {"enabled": False},
    "runcmd": {"remage": "remage"},
    "nersc": {"dvs_ro": False, "scratch": False},
}

#: Tier settings, ``simprod/config/tier/{tier}/{experiment}/settings``, keyed by
#: tier name.
TIER_SETTINGS: dict = {
    "opt": {
        "light_source": "optmap",
        "store_expected_pes": False,
        "scintillators": None,
    },
    "hit": {
        "simulate_psd": True,
        "simulate_psd_with_psl": False,
        "two_pass_aoe_correction": False,
        "tune_hpge_impurities_on_data": False,
        "store_surface_steps": False,
        "surface_depth_in_mm": 3,
    },
    "evt": {
        "skip_opt": False,
        "skip_hit": False,
        "add_random_coincidences": False,
        "random_coincidence_mode": "forced_trigger",
        "random_coincidence_runid": None,
    },
    "pdf": {
        "detector_groups": None,
        "fccd_tag_energy_thr_kev": None,
    },
}

#: HPGe parameter settings,
#: ``simprod/config/pars/{experiment}/geds/{par}/settings``, keyed by parameter
#: name.
PAR_SETTINGS: dict = {
    "aoemeanmod": {
        "simid_regex": "sis*_z*_slot*_Pb212_to_Pb208",
    },
    "elecmod": {
        "angle": "000",
        "sigma_start": 10.0,
        "tau_start": 50.0,
        "sigma_limits": (1.0, 50.0),
        "tau_limits": (1.0, 200.0),
        "comparison_window": (-500.0, 500.0),
        "plot_window": (-600.0, 600.0),
        "weight_power": 2.0,
        "max_calls": 1000,
        "dt_range_tuning": (600.0, 3000.0),
        "max_num_superpulses": 5,
        "minimiser_mode": "simplex",
        # None: taken from l200data, or from the built-in values of
        # utils.lookup_curr_mw_pars
        "curr_mwa": None,
    },
    "ssd": {
        "grid_size_in_mm": 0.5,
        "ssd_refinement_limits": [0.2, 0.1, 0.05, 0.02],
        "padding": 3,
    },
    "impurityscan": {
        # None: use grid_size_in_mm of the ssd settings
        "ssd": {"grid_size_in_mm": None},
        "drift_time_scan": {
            "min_energy_in_keV": 1500,
            "max_events": None,
        },
        "impurity_fit": {
            "energy_range_in_keV": [1500, 2500],
            "drift_time_uncertainty_in_ns": 50,
            "waveform_uncertainty": 0.5,
            "drift_time_observables": {
                "percentile": 90,
                "smoothing_in_ns": 50,
                "min_peak_fraction": 0.25,
            },
        },
    },
    "modeling": {
        "min_voltage_above_depletion_in_V": 100,
    },
    "superpulses": {
        "build_per_runid": False,
        "min_number_wfs": 10,
        "target_wfs": 100,
        "chi2_threshold": 3,
        "max_files": None,
        "tier_name": "auto",
        "charge_output": "wf_pz_win",
        "curr_output": "curr_av",
        "energy_output": "cuspEmax",
        "drift_time_slices": "1000:200:2000",
        "t0_field": "auto",
        "end_time_field": "geds/psd/low_aoe/time",
    },
}

#: Defaults of the ``elecmod`` settings in the impurity-curve scan
#: (``extract_elecmod_scan``), which reads the same settings file.
ELECMOD_SCAN_SETTINGS: dict = PAR_SETTINGS["elecmod"] | {
    "tau_start": 30.0,
    "tau_limits": (2.0, 100.0),
    "max_calls": 100,
}


def merge(user: Mapping | None, default: Mapping) -> AttrsDict:
    """Return a copy of `default` updated with the values in `user`.

    Nested blocks are merged recursively, any other value in `user` (lists
    included) replaces the default.
    """
    result = copy.deepcopy(dict(default))
    for k, v in (user or {}).items():
        if isinstance(result.get(k), Mapping) and isinstance(v, Mapping):
            result[k] = merge(v, result[k])
        else:
            result[k] = v
    return AttrsDict(result)
