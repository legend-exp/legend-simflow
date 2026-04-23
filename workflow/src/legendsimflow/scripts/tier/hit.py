# ruff: noqa: I002

# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>,
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

import argparse

import awkward as ak
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lgdo
import numpy as np
import pint
import pyg4ometry
import pygeomhpges
import pygeomtools
import reboost.hpge.psd
import reboost.hpge.surface
import reboost.hpge.utils
import reboost.math.functions
import reboost.math.stats
import reboost.spms
from dbetto import AttrsDict
from dbetto.utils import load_dict
from lgdo import lh5
from lgdo.lh5 import LH5Iterator
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import hpge_pars, nersc, patterns, utils
from legendsimflow import metadata as mutils
from legendsimflow import reboost as reboost_utils
from legendsimflow.metadata import get_tier_settings
from legendsimflow.profile import make_profiler
from legendsimflow.scripts import log_script_invocation
from legendsimflow.tcm import build_tcm


@snakemake_compatible(
    mapping={
        "stp_file": "input.stp_file",
        "jobid": "wildcards.jobid",
        "hit_file": "output[0]",
        "geom_file": "input.geom",
        "dtmap_files": "input.hpge_dtmaps",
        "currmod_files": "input.hpge_currmods",
        "simstat_part_file": "input.simstat_part_file",
        "detector_usabilities_file": "input.detector_usabilities[0]",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(description="Build the hit tier.")
    parser.add_argument("--stp-file", required=True, help="input stp tier file")
    parser.add_argument("--jobid", required=True, help="job ID wildcard")
    parser.add_argument("--hit-file", required=True, help="output hit tier file")
    parser.add_argument("--geom-file", required=True, help="GDML geometry file")
    parser.add_argument(
        "--dtmap-files", nargs="*", default=[], help="HPGe drift time map files"
    )
    parser.add_argument(
        "--currmod-files",
        nargs="*",
        default=[],
        help="HPGe current model parameter files",
    )
    parser.add_argument(
        "--simstat-part-file",
        required=True,
        help="simulation statistics partition file",
    )
    parser.add_argument(
        "--detector-usabilities-file",
        required=True,
        help="detector usabilities YAML file",
    )
    parser.add_argument("--log-file", default=None, help="log file")
    parser.add_argument(
        "--simflow-config",
        "--config",
        dest="simflow_config",
        required=True,
        help="simflow config YAML path",
    )
    args = parser.parse_args()

    config = utils.init_simflow_context(args.simflow_config, workflow=None).config

    stp_file = nersc.dvs_ro(config, args.stp_file)
    jobid = args.jobid
    hit_file = args.hit_file
    gdml_file = nersc.dvs_ro(config, args.geom_file)
    log_file = args.log_file
    metadata = config.metadata
    simstat_part_file = nersc.dvs_ro(config, args.simstat_part_file)
    l200data = config.paths.get("l200data", None)
    usabilities = AttrsDict(
        load_dict(nersc.dvs_ro(config, args.detector_usabilities_file))
    )

    # default resolutions/cuts for non-ON detectors, sourced from hit tier settings
    tier_hit_settings = get_tier_settings(config, "hit")
    dead_layer_fraction = tier_hit_settings.dead_layer_fraction
    buffer_len = tier_hit_settings.buffer_len
    eresmod_default = hpge_pars.build_energy_res_func_from_entry(
        tier_hit_settings.eresmod_default
    )
    aoeresmod_default = hpge_pars.build_aoe_res_func_from_entry(
        tier_hit_settings.aoeresmod_default
    )
    psdcuts_default = tier_hit_settings.psdcuts_default.to_dict()

    hit_file, move2cfs = nersc.make_on_scratch(config, hit_file)

    u = pint.UnitRegistry()

    # setup logging
    log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "tier-hit", parser, args)
    perf_block, print_perf, _ = make_profiler()

    # load the geometry and retrieve registered sensitive volume tables
    with perf_block("load_pygeom()"):
        geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
    sens_tables = pygeomtools.detectors.get_all_senstables(geom)

    # determine list of stp tables in the stp file
    ondisk_stp_tables = {}
    for tbl in lh5.ls(stp_file, "stp/*"):
        # ignore table names that start with underscore
        if not tbl.removeprefix("stp/").startswith("_"):
            ondisk_stp_tables[tbl] = False

    partitions = load_dict(simstat_part_file)[f"job_{jobid}"]

    # load TCM, to be used to chunk the event statistics according to the run partitioning
    log.debug("loading TCM")
    tcm = lh5.read_as("tcm", stp_file, library="ak")

    # extract the detector origin
    det_loc = lh5.read("detector_origins", stp_file)
    # convert to right format and add units
    # FIXME: units should be already present, to be fixed in remage
    det_loc = {
        k: [v[field].value for field in ("xloc", "yloc", "zloc")] * u.m
        for k, v in det_loc.items()
    }

    # loop over the partitions for this file
    for runid_idx, (runid, evt_idx_range) in enumerate(partitions.items()):
        log.info(
            "processing partition corresponding to %s [%d/%d], event range %s",
            runid,
            runid_idx + 1,
            len(partitions),
            evt_idx_range,
        )

        log.debug("loading energy resolution parameters")
        eresmod_pars_file = patterns.output_eresmod_filename(config, runid=runid)
        eresmod_pars_all = load_dict(eresmod_pars_file)
        energy_res_func = hpge_pars.build_energy_res_func_dict(
            l200data,
            metadata,
            runid,
            energy_res_pars=eresmod_pars_all,
        )  # FWHM

        log.debug("loading A/E resolution parameters")
        aoeresmod_pars_file = patterns.output_aoeresmod_filename(config, runid=runid)
        aoeresmod_pars_all = load_dict(aoeresmod_pars_file)
        aoe_res_func = hpge_pars.build_aoe_res_func_dict(
            l200data,
            metadata,
            runid,
            aoe_res_pars=aoeresmod_pars_all,
        )

        log.debug("loading current pulse model parameters")
        currmod_pars_file = patterns.output_currmod_merged_filename(config, runid=runid)
        currmod_pars_all = AttrsDict(load_dict(currmod_pars_file))

        log.debug("loading PSD cut values")
        psdcuts_file = patterns.output_psdcuts_filename(config, runid=runid)
        psdcuts_all = load_dict(psdcuts_file)

        # loop over the sensitive volume tables registered in the geometry
        for det_idx, (det_name, geom_meta) in enumerate(sens_tables.items()):
            log.debug(
                "looking for data from sensitive volume table %s (uid=%s)...",
                det_name,
                geom_meta.uid,
            )

            stp_table_name = f"stp/{det_name}"

            # only process the HPGe output
            if geom_meta.detector_type != "germanium":
                # no bookkeeping of this table
                ondisk_stp_tables.pop(stp_table_name, None)
                continue

            if stp_table_name not in ondisk_stp_tables:
                log.warning(
                    "detector %s not found in %s. "
                    "possibly because it was not read-out or there were no hits recorded",
                    det_name,
                    stp_file,
                )
                continue

            # get the usability
            det_info = usabilities[runid][det_name]
            if det_info is None:
                log.warning(
                    "usability not found for %s in %s, defaulting to on",
                    det_name,
                    runid,
                )
                usability = "on"
                psd_usability = "valid"
            else:
                usability = det_info.usability
                psd_usability = det_info.psd_usability
            psd_usability_code = mutils.encode_psd_usability(psd_usability)

            log.debug("looking for indices of hit table rows to read...")
            with perf_block("get_remage_hit_range()"):
                i_start, n_entries = reboost_utils.get_remage_hit_range(
                    tcm, det_name, geom_meta.uid, evt_idx_range
                )

            # initialize the stp file iterator
            # NOTE: if the entry list is empty, there will be no processing but an
            # empty output table will be nonetheless created. this is important for
            # the buil_tcm() step at the end

            iterator = LH5Iterator(
                stp_file,
                stp_table_name,
                i_start=i_start,
                n_entries=n_entries,
                buffer_len=buffer_len,
            )

            log.info(
                "processing the %s output table [%d/%d]...",
                det_name,
                det_idx + 1,
                len(sens_tables),
            )

            log.debug("creating an pygeomhpges.HPGe object")
            pyobj = pygeomhpges.make_hpge(
                geom_meta.metadata, registry=None, allow_cylindrical_asymmetry=False
            )

            fccd = mutils.get_sanitized_fccd(metadata, det_name)

            # NOTE: we don't use the script arg but we use the (known) file patterns. more robust
            dt_map = reboost_utils.load_hpge_dtmaps(config, det_name, runid)

            # load parameters of the current model
            pars = currmod_pars_all.get(det_name, None)
            currmod_pars = (
                pars.get("current_pulse_pars", None) if pars is not None else None
            )

            can_model_psd = dt_map is not None and currmod_pars is not None

            if not can_model_psd and usability == "on" and psd_usability == "valid":
                log.warning(
                    "%s is ON with valid PSD in data but its PSD response "
                    "could not be simulated (drift-time map or current model missing).",
                    det_name,
                )

            # iterate over input data
            for lgdo_chunk in iterator:
                chunk = lgdo_chunk.view_as("ak", with_units=True)

                with perf_block("distance_to_surface()"):
                    _distance_to_nplus = reboost.hpge.surface.distance_to_surface(
                        chunk.xloc,
                        chunk.yloc,
                        chunk.zloc,
                        pyobj,
                        det_loc[det_name],
                        distances_precompute=chunk.dist_to_surf,
                        precompute_cutoff=(fccd + 1),
                        surface_type="nplus",
                    )

                with perf_block("piecewise_linear_activeness()"):
                    _activeness = reboost.math.functions.piecewise_linear_activeness(
                        _distance_to_nplus,
                        fccd_in_mm=fccd,
                        dlf=dead_layer_fraction,
                    )

                edep_active = chunk.edep * _activeness
                energy_true = ak.sum(edep_active, axis=-1)

                # Validation counterpart to the collection in
                # extract_hpge_observables_models: ON detectors must have curves
                # (hard error); others fall back to eresmod_default (soft warning).
                # TODO: move to a separate function to clean up

                if det_name in energy_res_func:
                    energy_res = energy_res_func[det_name](energy_true)

                elif usability == "on":
                    msg = f"{det_name} is ON but no energy resolution curves are available"
                    raise RuntimeError(msg)
                else:
                    log.warning(
                        "%s is marked as '%s' and no energy resolution curves are available. "
                        "using eresmod_default from hit tier settings",
                        det_name,
                        usability,
                    )
                    energy_res = eresmod_default(energy_true)

                if det_name in aoe_res_func:
                    aoe_res = aoe_res_func[det_name](energy_true)
                elif usability == "on" and psd_usability != "missing" and can_model_psd:
                    msg = (
                        f"{det_name} is ON with valid PSD but no A/E resolution "
                        "curves are available"
                    )
                    raise RuntimeError(msg)
                else:
                    log.warning(
                        "%s is marked as '%s' (psd_usability='%s') "
                        "and no A/E resolution curves are available. using aoeresmod_default",
                        det_name,
                        usability,
                        psd_usability,
                    )
                    aoe_res = aoeresmod_default(energy_true)

                if det_name in psdcuts_all:
                    psdcuts = AttrsDict(
                        utils.sanitize_dict_with_defaults(
                            psdcuts_all[det_name],
                            psdcuts_default,
                        )
                    )
                elif usability == "on" and psd_usability != "missing" and can_model_psd:
                    msg = f"{det_name} is ON with valid PSD but no PSD cut values are available"
                    raise RuntimeError(msg)
                else:
                    log.warning(
                        "%s is marked as '%s' (psd_usability='%s') "
                        "and no PSD cut values are available. using psdcuts_default",
                        det_name,
                        usability,
                        psd_usability,
                    )
                    psdcuts = AttrsDict(psdcuts_default)

                # smear energy with detector resolution
                energy = reboost_utils.gauss_smear(
                    energy_true,
                    energy_res / 2.35482,
                )

                # PSD: if the drift time map is none, it means that we don't
                # have the detector model to simulate PSD in a more advanced
                # way

                # default to NaN
                _drift_time = ak.full_like(chunk.xloc, fill_value=np.nan)
                aoe = np.full(len(chunk), np.nan)
                aoe_class = np.full(len(chunk), np.nan)
                is_single_site = np.full(len(chunk), False)
                t_max = np.full(len(chunk), np.nan)

                if dt_map is not None and currmod_pars is not None:
                    log.info("computing PSD observables")

                    with perf_block("hpge_corrected_drift_time()"):
                        _drift_time = reboost_utils.hpge_corrected_drift_time(
                            chunk, dt_map, det_loc[det_name]
                        )
                    utils.check_nans_leq(
                        _drift_time, "_drift_time", 0.01, min_entries=1000
                    )

                    with perf_block("hpge_max_current()"):
                        _a_max_true = reboost_utils.hpge_max_current(
                            edep_active, _drift_time, currmod_pars
                        )
                    utils.check_nans_leq(
                        _a_max_true, "_a_max_true", 0.01, min_entries=1000
                    )

                    # Apply current resolution smearing based on configured A/E noise parameters
                    _a_max = reboost_utils.gauss_smear(
                        _a_max_true, pars.current_reso / pars.mean_aoe
                    )

                    # gauss_smear floors negative-smeared energy to np.finfo(float).tiny;
                    # dividing by that sentinel produces values that overflow float32.
                    has_physical_energy = energy > np.finfo(np.float64).tiny
                    aoe = ak.where(has_physical_energy, _a_max / energy, np.nan)

                    # ...and A/E classifier
                    # NOTE: we use the resolution determined from data here instead
                    # of the intrinsic simulated ones due to noise
                    aoe_class = ak.where(
                        has_physical_energy, (aoe - 1) / aoe_res, np.nan
                    )

                    # ...and PSD flag
                    is_single_site = (aoe_class > psdcuts.aoe.low_side) & (
                        aoe_class < psdcuts.aoe.high_side
                    )

                    # also calculate drift time at A position
                    # FIXME: this is wasting compute resources, max_current should
                    # return (maxA, t_maxA)
                    with perf_block("hpge_max_current()"):
                        t_max = reboost_utils.hpge_max_current(
                            edep_active,
                            _drift_time,
                            currmod_pars,
                            return_mode="max_time",
                        )

                out_table = reboost_utils.make_output_chunk(lgdo_chunk)

                out_table.add_field(
                    "energy",
                    lgdo.Array(
                        np.asarray(energy, dtype=np.float32), attrs={"units": "keV"}
                    ),
                )
                out_table.add_field(
                    "drift_time_amax",
                    lgdo.Array(
                        np.asarray(t_max, dtype=np.float32), attrs={"units": "ns"}
                    ),
                )
                out_table.add_field(
                    "aoe_raw", lgdo.Array(np.asarray(aoe, dtype=np.float32))
                )
                out_table.add_field(
                    "aoe", lgdo.Array(np.asarray(aoe_class, dtype=np.float32))
                )
                out_table.add_field("is_single_site", lgdo.Array(is_single_site))

                _, period, run, _ = mutils.parse_runid(runid)
                field_vals = [period, run, mutils.encode_usability(usability)]
                for i, field in enumerate(["period", "run", "usability"]):
                    out_table.add_field(
                        field,
                        lgdo.Array(np.full(shape=len(chunk), fill_value=field_vals[i])),
                    )
                out_table.add_field(
                    "psd_usability",
                    lgdo.Array(
                        np.full(shape=len(chunk), fill_value=psd_usability_code)
                    ),
                )

                with perf_block("write_chunk()"):
                    reboost_utils.write_chunk(
                        out_table,
                        f"/hit/{det_name}",
                        hit_file,
                        geom_meta.uid,
                    )

            # this table has been processed
            ondisk_stp_tables[stp_table_name] = True

    # sanity check that all the stp tables were processed. this is important for
    # TCM consistency across all tiers later on
    not_done = [k for k, v in ondisk_stp_tables.items() if not v]
    if not_done:
        msg = f"stp tables {not_done} were not processed!"
        raise RuntimeError(msg)

    log.debug("building the TCM")
    build_tcm(hit_file, hit_file)

    with perf_block("move_to_cfs()"):
        move2cfs()

    print_perf()


if __name__ == "__main__":
    main()
