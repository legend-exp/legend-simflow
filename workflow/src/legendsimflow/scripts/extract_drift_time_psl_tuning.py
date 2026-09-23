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
import logging

import awkward as ak
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import lh5
import pint
import pyg4ometry
import pygeomhpges
import pygeomtools
import reboost.hpge
import reboost.hpge.surface
import reboost.hpge.utils
import reboost.math
from lgdo import Array, Scalar, Struct
from lh5 import LH5Iterator
from lh5.io import read_n_rows
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import metadata as mutils
from legendsimflow import nersc, psl, utils
from legendsimflow.hpge_electronics_tuning import load_elecmod
from legendsimflow.psl import validate_ssd_scan_grid
from legendsimflow.reboost import cluster_steps, get_rz, mask_with_units
from legendsimflow.scripts import log_script_invocation


@snakemake_compatible(
    mapping={
        "stp_files": "input.stp_files",
        "hpge_detector": "wildcards.hpge_detector",
        "drift_time_file": "output",
        "elecmod": "input.elecmod",
        "geom_file": "input.geom",
        "psl_file": "input.psl_file",
        "log_file": "log[0]",
        "simflow_config": "config",
        "max_events": "params.max_events",
        "energy_cut": "params.energy_cut",
    }
)
def main() -> None:
    """Compute event drift times over a scan of pulse-shape libraries.

    The scan covers a two-dimensional grid of impurity-curve slope and
    depletion voltage, in the format described in
    :func:`legendsimflow.psl.validate_ssd_scan_grid`. The steps are:

    1. Read the geometry of the HPGe detector from the GDML file
       (``--geom-file``).
    2. Read the ideal pulse-shape libraries of the scan from ``--psl-file``
       (:func:`legendsimflow.psl.load_ideal_psl_scan`) and convolve each of
       them with the electronics response built from the parameters in
       ``--elecmod`` (:func:`legendsimflow.psl.convolve_elecmod_scan`).
    3. Read the steps from ``--stp-files`` in chunks, cluster them, and keep
       the events depositing more than 1500 keV in the active volume.
    4. For each grid point, interpolate the drift time of every cluster
       (:func:`reboost.hpge.drift_time_crystal_axes`) and take the drift time
       of the event as the one of maximum current
       (:func:`reboost.hpge.maximum_current`).

    Processing stops once ``--max-events`` events have passed the energy cut.

    The output file ``--drift-time-file`` follows the same scan format as the
    input, with a ``drift_time`` array in ns at each grid point. Next to
    ``psl_scan`` it carries an ``energy`` array with the energy each event
    deposits in the active volume, in keV, in the same order as the drift
    times.

    The file also contains the ``grid_info`` structure with the grid definition, and a ``weight``
    scalar with the fraction of the total available MC statistics that were used for the output files.
    """
    parser = argparse.ArgumentParser(
        description="Compute event drift times over a scan of pulse-shape libraries."
    )
    parser.add_argument(
        "--stp-files", nargs="+", required=True, help="Input stp files."
    )
    parser.add_argument(
        "--drift-time-file", required=True, help="output drift time file"
    )
    parser.add_argument("--hpge-detector", required=True, help="HPGe detector name")
    parser.add_argument(
        "--psl-file",
        required=True,
        help="scan of ideal HPGe pulse-shape libraries.",
    )
    parser.add_argument("--elecmod", required=True, help="HPGe electronics model file.")

    parser.add_argument("--geom-file", required=True, help="input geom file")
    parser.add_argument("--simflow-config", help="simflow config file")
    parser.add_argument("--log-file", help="log file")
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="stop after this many events pass the energy cut (default: use all)",
    )

    parser.add_argument(
        "--energy-cut",
        type=float,
        default=1500.0,
        help="energy cut in keV (default: 1500 keV)",
    )
    args = parser.parse_args()
    det = args.hpge_detector
    log_file = args.log_file
    files = args.stp_files

    # get file paths
    if args.simflow_config is not None:
        config = utils.init_simflow_context(args.simflow_config, workflow=None).config
        dt_file, move2cfs = nersc.make_on_scratch(config, args.drift_time_file)
        metadata = config.metadata

        log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
        log_script_invocation(log, "extract-drift-time-psl-tuning", parser, args)

        files = [nersc.dvs_ro(config, s) for s in files]
        gdml_file = nersc.dvs_ro(config, args.geom_file)
        psl_file = nersc.dvs_ro(config, args.psl_file)
        elecmod_file = nersc.dvs_ro(config, args.elecmod)

    else:
        dt_file = args.drift_time_file

        def move2cfs():
            return None

        metadata = None

        logging.basicConfig(
            level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
        )
        log = logging.getLogger(__name__)

        gdml_file = args.geom_file
        psl_file = args.psl_file
        elecmod_file = args.elecmod

    # other setup
    u = pint.get_application_registry()

    # setup logging
    perf_block, print_perf, _ = reboost.make_profiler()

    # get the geometry
    with perf_block("load_pygeom()"):
        geom = pyg4ometry.gdml.Reader(gdml_file).getRegistry()
        sens_tables = pygeomtools.detectors.get_all_senstables(geom)

    # get the files
    det_loc = lh5.read("detector_origins", files[0])
    det_loc = {
        k: [v[field].value for field in ("xloc", "yloc", "zloc")] * u.m
        for k, v in det_loc.items()
    }

    stp_table_name = f"stp/{det}"
    geom_meta = sens_tables[det]

    buffer_len = 200000
    if args.max_events is not None:
        buffer_len = min(args.max_events, buffer_len)

    iterator = LH5Iterator(
        files,
        stp_table_name,
        i_start=0,
        buffer_len=buffer_len,
    )

    # extract necessary geometry information
    pyobj = pygeomhpges.make_hpge(
        geom_meta.metadata, registry=None, allow_cylindrical_asymmetry=False
    )

    if metadata is not None:
        fccd = mutils.get_sanitized_fccd(metadata, det)
    else:
        msg = (
            f"No metadata found in config file. Using default FCCD of 1.0 mm for {det}."
        )
        log.info(msg)
        fccd = 1.0

    log.info("... extracting stps and energies")
    # loop over steps
    n_read = 0
    n_used = 0
    stps = []

    for lgdo_chunk in iterator:
        chunk = lgdo_chunk.view_as("ak", with_units=True)

        # remove events with energy below ECUT
        n_read += len(chunk)
        chunk = mask_with_units(chunk, ak.sum(chunk.edep, axis=-1) > args.energy_cut)

        log.info("... cluster steps")
        # cluster steps
        with perf_block("cluster_steps()"):
            chunk_new = cluster_steps(
                chunk, surf_cut=2, threshold_in_mm=1, threshold_surf_in_mm=0.05
            )
        log.info("... compute energy")

        # add some clustering
        with perf_block("activeness"):
            _distance_to_nplus = reboost.hpge.surface.distance_to_surface(
                chunk_new.xloc,
                chunk_new.yloc,
                chunk_new.zloc,
                pyobj,
                det_loc[det],
                distances_precompute=chunk_new.dist_to_surf,
                precompute_cutoff=(fccd + 1),
                surface_type="nplus",
            )

            _activeness = reboost.math.piecewise_linear_activeness(
                _distance_to_nplus,
                fccd_in_mm=fccd,
                dlf=0.5,
            )

            edep_active = chunk_new.edep * _activeness
            energy_true = ak.sum(edep_active, axis=-1)

            # cut sub threshold events
            chunk_new = mask_with_units(chunk_new, energy_true > args.energy_cut)
            edep_active = edep_active[energy_true > args.energy_cut]
            energy_true = energy_true[energy_true > args.energy_cut]

        # now get drift times
        log.info("... get drift times")

        chunk_new["energy_true"] = energy_true
        chunk_new["edep_active"] = edep_active

        stps.append(chunk_new)
        n_used += len(energy_true)
        if args.max_events is not None and n_used >= args.max_events:
            log.info(
                "... we have now gathered %d events and we only needed %d so we are stopping",
                n_used,
                args.max_events,
            )
            break

    stps = ak.concatenate(stps)

    lh5.write(
        Struct({"energy", Array(stps.energy_true)}),
        f"/{det}/energy",
        dt_file,
        wo_mode="of",
    )

    log.info("... start extraction of drift times for %d events", n_used)

    sigma, tau = load_elecmod(elecmod_file)

    # lookup the psl scan groups and grid info
    psl_groups, grid_info = psl.lookup_ideal_psl_scan_groups(psl_file)

    for slope, depv_psls in psl_groups.items():
        for depv in depv_psls:
            with perf_block("read_psl"):
                ideal_psl = lh5.read(depv_psls[depv], psl_file)

                realistic_psl, dt_maps = psl.convolve_elecmod(
                    ideal_psl, sigma=sigma, tau=tau
                )

            with perf_block("drift_time"):
                _drift_time = reboost.hpge.drift_time_crystal_axes(
                    stps.xloc,
                    stps.yloc,
                    stps.zloc,
                    dt_maps,
                    coord_offset=det_loc[det],
                )
                _r, _z = get_rz(det_loc[det], stps)

                drift_time = {
                    "drift_time": reboost.hpge.maximum_current(
                        stps.edep_active,
                        _drift_time,
                        times=None,
                        r=_r,
                        z=_z,
                        template=realistic_psl,
                        return_mode="max_time",
                    )
                }
                out = Array(drift_time)
                lh5.write(
                    out,
                    f"/{det}/psl_scan/{slope}/{depv}",
                    dt_file,
                    wo_mode="append_column",
                )

    n_tot = sum([read_n_rows(f"stp/{det}", stp_file) for stp_file in files])
    log.info(
        "computed drift times for %d events above %d keV, out of %d read out of a possible %d events",
        n_used,
        args.energy_cut,
        n_read,
        n_tot,
    )
    # save a weight (used for normalisation)
    weight = n_read / n_tot
    lh5.write(Scalar(weight), f"{det}/weight", dt_file, wo_mode="append_column")

    # write the grid definition
    lh5.write(grid_info, f"{det}/grid_info/", dt_file, wo_mode="append_column")

    is_valid = validate_ssd_scan_grid(dt_file, det)
    if not is_valid:
        msg = f"output file {dt_file} does not have a valid scan grid"
        raise ValueError(msg)

    with perf_block("move_to_cfs()"):
        move2cfs()

    print_perf()


if __name__ == "__main__":
    main()
