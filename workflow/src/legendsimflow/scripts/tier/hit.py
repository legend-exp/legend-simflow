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

import awkward as ak
import dbetto.utils
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
from lgdo import lh5
from lgdo.lh5 import LH5Iterator

from legendsimflow import hpge_pars, nersc, patterns, utils
from legendsimflow import metadata as mutils
from legendsimflow import reboost as reboost_utils

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

stp_file = args.input.stp_file
jobid = args.wildcards.jobid
hit_file = args.output[0]
gdml_file = args.input.geom
log_file = args.log[0]
metadata = args.config.metadata
hpge_dtmap_files = args.input.hpge_dtmaps
hpge_currmods_files = args.input.hpge_currmods
hpge_eresmods_files = args.input.hpge_eresmods
simstat_part_file = args.input.simstat_part_file
l200data = args.config.paths.l200data

BUFFER_LEN = "500*MB"

u = pint.UnitRegistry()


def DEFAULT_ENERGY_RES_FUNC(energy):
    return 2.5 * np.sqrt(energy / 2039)  # FWHM


# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)

# load the geometry and retrieve registered sensitive volume tables
geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
sens_tables = pygeomtools.detectors.get_all_senstables(geom)

# determine list of stp tables in the stp file
ondisk_stp_tables = {}
for tbl in lh5.ls(stp_file, "stp/*"):
    # ignore table names that start with underscore
    if not tbl.removeprefix("stp/").startswith("_"):
        ondisk_stp_tables[tbl] = False

partitions = dbetto.utils.load_dict(simstat_part_file)[f"job_{jobid}"]

# load TCM, to be used to chunk the event statistics according to the run partitioning
msg = "loading TCM"
log.debug(msg)
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
    msg = (
        f"processing partition corresponding to {runid} "
        f"[{runid_idx + 1}/{len(partitions)}], event range {evt_idx_range}"
    )
    log.info(msg)

    msg = "loading energy resolution parameters"
    log.debug(msg)
    eresmod_pars_file = patterns.output_eresmod_filename(
        snakemake.config,  # noqa: F821
        runid=runid,
    )
    eresmod_pars_all = dbetto.utils.load_dict(eresmod_pars_file)
    energy_res_func = hpge_pars.build_energy_res_func_dict(
        l200data,
        metadata,
        runid,
        energy_res_pars=eresmod_pars_all,
    )  # FWHM

    msg = "loading current pulse model parameters"
    log.debug(msg)
    currmod_pars_file = patterns.output_currmod_merged_filename(
        snakemake.config,  # noqa: F821
        runid=runid,
    )
    currmod_pars_all = AttrsDict(dbetto.utils.load_dict(currmod_pars_file))

    # loop over the sensitive volume tables registered in the geometry
    for det_idx, (det_name, geom_meta) in enumerate(sens_tables.items()):
        msg = f"looking for data from sensitive volume table {det_name} (uid={geom_meta.uid})..."
        log.debug(msg)

        stp_table_name = f"stp/{det_name}"

        # only process the HPGe output
        if geom_meta.detector_type != "germanium":
            # no bookkeeping of this table
            ondisk_stp_tables.pop(stp_table_name, None)
            continue

        if stp_table_name not in ondisk_stp_tables:
            msg = (
                f"detector {det_name} not found in {stp_file}. "
                "possibly because it was not read-out or there were no hits recorded"
            )
            log.warning(msg)

            continue

        # get the usability
        usability = mutils.usability(metadata, det_name, runid=runid, default="on")

        msg = "looking for indices of hit table rows to read..."
        log.debug(msg)
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
            buffer_len=BUFFER_LEN,
        )

        msg = f"processing the {det_name} output table [{det_idx + 1}/{len(sens_tables)}]..."
        log.info(msg)

        log.debug("creating an pygeomhpges.HPGe object")
        pyobj = pygeomhpges.make_hpge(
            geom_meta.metadata, registry=None, allow_cylindrical_asymmetry=False
        )

        fccd = mutils.get_sanitized_fccd(metadata, det_name)

        # NOTE: we don't use the script arg but we use the (known) file patterns. more robust
        dt_map = reboost_utils.load_hpge_dtmaps(snakemake.config, det_name, runid)  # noqa: F821

        # load parameters of the current model
        pars = currmod_pars_all.get(det_name, None)
        currmod_pars = (
            pars.get("current_pulse_pars", None) if pars is not None else None
        )

        n_tot = 0
        # iterate over input data
        for lgdo_chunk in iterator:
            chunk = lgdo_chunk.view_as("ak", with_units=True)

            n_tot += len(chunk)

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

            _activeness = reboost.math.functions.piecewise_linear_activeness(
                _distance_to_nplus,
                fccd_in_mm=fccd,
                dlf=0.5,
            )

            edep_active = chunk.edep * _activeness
            energy_true = ak.sum(edep_active, axis=-1)

            # smear energy with detector resolution
            if det_name in energy_res_func:
                energy_res = energy_res_func[det_name](energy_true)
            elif usability != "off":
                msg = (
                    f"{det_name} is marked as '{usability}' but no "
                    "resolution curves are available. this is unacceptable!"
                )
                raise RuntimeError(msg)
            else:
                msg = (
                    f"{det_name} is marked as '{usability}' but no "
                    "resolution curves are available. using default "
                    "resolution of 2.5 keV FWHM at 2 MeV!"
                )
                log.warning(msg)
                energy_res = DEFAULT_ENERGY_RES_FUNC(energy_true)

            energy = reboost.math.stats.gaussian_sample(
                energy_true,
                energy_res / 2.35482,
            )

            # energy can't be negative as a result of smearing and later we
            # divide for it
            energy = ak.where(
                (energy <= 0) & (energy_true >= 0), np.finfo(float).tiny, energy
            )

            # PSD: if the drift time map is none, it means that we don't
            # have the detector model to simulate PSD in a more advanced
            # way

            # default to NaN
            drift_time = ak.full_like(chunk.xloc, fill_value=np.nan)
            aoe = np.full(len(chunk), np.nan)

            if dt_map is not None and currmod_pars is not None:
                msg = "computing PSD observables"
                log.info(msg)

                drift_time = reboost_utils.hpge_corrected_drift_time(
                    chunk, dt_map, det_loc[det_name]
                )
                # TODO: fix dtmap nan issue
                utils.check_nans_leq(drift_time, "drift_time", 0.9)

                _a_max = reboost_utils.hpge_max_current(
                    edep_active, drift_time, currmod_pars
                )
                # Apply current resolution smearing based on configured A/E noise parameters
                a_sigma = pars.current_reso / pars.mean_aoe

                _a_max = reboost.math.stats.gaussian_sample(_a_max, sigma=a_sigma)
                # TODO: fix dtmap nan issue
                utils.check_nans_leq(_a_max, "max_current", 0.9)

                aoe = _a_max / energy

            out_table = reboost_utils.make_output_chunk(lgdo_chunk)

            out_table.add_field("energy", lgdo.Array(energy, attrs={"units": "keV"}))
            out_table.add_field(
                "drift_time", lgdo.VectorOfVectors(drift_time, attrs={"units": "ns"})
            )
            out_table.add_field("aoe", lgdo.Array(aoe))

            _, period, run, _ = mutils.parse_runid(runid)
            field_vals = [period, run, mutils.encode_usability(usability)]
            for i, field in enumerate(["period", "run", "usability"]):
                out_table.add_field(
                    field,
                    lgdo.Array(np.full(shape=len(chunk), fill_value=field_vals[i])),
                )

            reboost_utils.write_chunk(
                out_table,
                f"/hit/{det_name}",
                hit_file,
                geom_meta.uid,
            )

        assert n_tot == n_entries
        # this table has been processed
        ondisk_stp_tables[stp_table_name] = True

# sanity check that all the stp tables were processed. this is important for
# TCM consistency across all tiers later on
not_done = [k for k, v in ondisk_stp_tables.items() if not v]
if not_done:
    msg = f"stp tables {not_done} were not processed!"
    raise RuntimeError(msg)

log.debug("building the TCM")
reboost_utils.build_tcm(hit_file, hit_file)
