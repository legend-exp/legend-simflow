# ruff: noqa: I002

# Copyright (C) 2026 Toby Dixon <toby.dixon.23@ucl.ac.uk>,
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
from pathlib import Path

import awkward as ak
import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
import numpy as np
import pint
import pyg4ometry
import pygeomhpges
import pygeomtools
import reboost.hpge.surface
import reboost.hpge.utils
from lgdo import Table
from lh5 import LH5Iterator
from reboost import units
from reboost.shape.cluster import apply_cluster, cluster_by_step_length
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import metadata as mutils
from legendsimflow import nersc, psl, utils
from legendsimflow import reboost as reboost_utils
from legendsimflow.scripts import log_script_invocation

PADDING = (0, 2000)
N_MAX = 10000


def mask_with_units(data: ak.Array, mask: ak.Array) -> ak.Array:
    """Mask an awkward array with units attached, preserving the units."""
    u = {field: units.get_unit_str(data[field]) for field in data.fields}
    data = data[mask]

    for field in data.fields:
        data[field] = units.attach_units(data[field], u[field])

    return data


def cluster_steps(chunk: ak.Array, **kwargs) -> ak.Array:
    """Cluster steps in a chunk of events, returning a new chunk with clustered steps."""
    clusters = cluster_by_step_length(
        ak.ones_like(chunk.trackid),
        chunk.xloc,
        chunk.yloc,
        chunk.zloc,
        units.units_conv_ak(chunk.dist_to_surf, "mm"),
        **kwargs,
    )

    xc = _apply_cluster(clusters, chunk.xloc, mode="mean")
    yc = _apply_cluster(clusters, chunk.yloc, mode="mean")
    zc = _apply_cluster(clusters, chunk.zloc, mode="mean")
    ec = _apply_cluster(clusters, chunk.edep, mode="sum")
    dc = _apply_cluster(clusters, chunk.dist_to_surf, mode="mean")

    return ak.Array(
        {"xloc": xc, "yloc": yc, "zloc": zc, "dist_to_surf": dc, "edep": ec}
    )


def _apply_cluster(clusters: ak.Array, data: ak.Array, mode: str = "sum") -> ak.Array:
    """Apply clustering to a data array, returning the clustered data.

    Parameters
    ----------
    clusters
        The cluster indices for each step in the data array.
    data
        The data array to be clustered.
    mode
        The mode of clustering to apply. Can be "sum" or "mean". Defaults to "sum".
    """
    unit = units.get_unit_str(data)

    data_cluster = apply_cluster(clusters, data)
    if mode == "sum":
        return units.attach_units(ak.sum(data_cluster, axis=-1), unit)
    if mode == "mean":
        return units.attach_units(ak.mean(data_cluster, axis=-1), unit)
    msg = f"Mode {mode} not recognised. Must be 'sum' or 'mean'."
    raise ValueError(msg)

    return units.attach_units(data_cluster, unit)


def get_rz(det_loc, chunk: ak.Array) -> tuple[ak.Array, ak.Array]:
    """Extract the r and z coordinates of a chunk of events, given the detector location."""
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

    return _r, _z


ECUT = 1500
N_MAX = 10000


@snakemake_compatible(
    mapping={
        "stp_file": "input.simid",
        "jobid": "wildcards.jobid",
        "hpge_detector": "wildcards.hpge_detector",
        "drift_time_file": "output",
        "elecmod": "input.elecmod",
        "geom_file": "input.geom",
        "psl_file": "input.psl_file",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(description="Build the hit tier.")
    parser.add_argument("--stp-file", required=True, help="input stp tier file")
    parser.add_argument(
        "--drift-time-file", required=True, help="output drift time file"
    )
    parser.add_argument("--hpge-detector", required=True, help="HPGe detector name")
    parser.add_argument(
        "--psl-file", required=True, help="HPGe realistic pulse shape library file."
    )
    parser.add_argument("--elecmod", required=True, help="HPGe electronics model file.")

    parser.add_argument("--det", required=True, help="detector")
    parser.add_argument("--geom-file", required=True, help="input geom file")
    parser.add_argument("--simflow-config", required=True, help="simflow config file")
    parser.add_argument("--log-file", required=True, help="log file")

    args = parser.parse_args()
    det = args.hpge_detector

    # get file paths
    config = utils.init_simflow_context(args.simflow_config, workflow=None).config
    stp_file = nersc.dvs_ro(config, args.stp_file)
    gdml_file = nersc.dvs_ro(config, args.geom_file)
    log_file = args.log_file

    dt_file, move2cfs = nersc.make_on_scratch(config, args.drift_time_file)

    # other setup
    metadata = config.metadata
    u = pint.UnitRegistry()

    # setup logging
    log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
    log_script_invocation(log, "extract-drift-time-psl-tuning", parser, args)
    perf_block, print_perf, _ = reboost.make_profiler()

    # get the geometry
    with perf_block("load_pygeom()"):
        geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
        sens_tables = pygeomtools.detectors.get_all_senstables(geom)

    # get the files
    files = Path(stp_file).glob("*.lh5")
    det_loc = lh5.read("detector_origins", files[0])
    det_loc = {
        k: [v[field].value for field in ("xloc", "yloc", "zloc")] * u.m
        for k, v in det_loc.items()
    }

    stp_table_name = f"stp/{det}"
    geom_meta = sens_tables[det]

    iterator = LH5Iterator(
        files,
        stp_table_name,
        i_start=0,
        buffer_len=200000,
    )

    # extract necessary geometry information
    pyobj = pygeomhpges.make_hpge(
        geom_meta.metadata, registry=None, allow_cylindrical_asymmetry=False
    )

    fccd = mutils.get_sanitized_fccd(metadata, det)

    with perf_block("load_psl()"):
        ideal_psls = psl.load_ideal_psl_scan(args.psl_file, det)
        elecmod_pars = dbetto.utils.load_dict(args.elecmod, det)["best_fit"]
        realistic_psl, psl_dt_maps = psl.convolve_elecmod_scan(
            ideal_psls, elecmod_pars, padding=PADDING
        )

    # loop over steps
    for lgdo_chunk in iterator:
        chunk = lgdo_chunk.view_as("ak", with_units=True)

        # remove events with energy below ECUT
        chunk = mask_with_units(chunk, ak.sum(chunk.edep, axis=-1) > ECUT)

        # cluster steps
        with perf_block("cluster_steps()"):
            chunk_new = cluster_steps(
                chunk, surf_cut=2, threshold_in_mm=1, threshold_surf_in_mm=0.05
            )

        # add some clustering
        with perf_block("activeness"):
            _distance_to_nplus = reboost.hpge.surface.distance_to_surface(
                chunk_new.xloc,
                chunk_new.xloc,
                chunk_new.xloc,
                pyobj,
                det_loc[det],
                distances_precompute=chunk_new.dist_to_surf,
                precompute_cutoff=(fccd + 1),
                surface_type="nplus",
            )

            _activeness = reboost.math.functions.piecewise_linear_activeness(
                _distance_to_nplus,
                fccd_in_mm=fccd,
                dlf=0.5,
            )

            edep_active = chunk_new.edep * _activeness
            energy_true = ak.sum(edep_active, axis=-1)

        # now get drift times

        with perf_block("drift_time"):
            drift_times = {}

            for idx, dep in enumerate(psl_dt_maps):
                if idx % 10 == 0:
                    msg = f"Processing drift time for {dep} ({idx}/{len(psl_dt_maps)})"
                    log.info(msg)

                _drift_time = reboost_utils.hpge_corrected_drift_time(
                    chunk_new, psl_dt_maps[dep], det_loc[det]
                )
                _r, _z = get_rz(det_loc[det], chunk_new)

                drift_times[dep] = reboost.hpge.psd.maximum_current(
                    edep_active,
                    _drift_time,
                    times=None,
                    r=_r,
                    z=_z,
                    template=realistic_psl[dep],
                    return_mode="max_time",
                )
        if drift_times == {}:
            out = Table(ak.Array({"energy": energy_true}))
        else:
            out = Table(ak.Array({"energy": energy_true, "psd": drift_times}))

        reboost_utils.write_chunk(
            out,
            f"/drift_time/{det}",
            dt_file,
            geom_meta.uid,
        )

    with perf_block("move_to_cfs()"):
        move2cfs()

    print_perf()


if __name__ == "__main__":
    main()
