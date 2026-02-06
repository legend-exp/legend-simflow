# Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
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

import inspect
import os
from datetime import datetime
from pathlib import Path

import legenddataflowscripts as ldfs
import pyg4ometry as pg4
from dbetto import AttrsDict
from legendmeta import LegendMetadata
from lgdo import lh5
from numpy.typing import ArrayLike
from reboost.hpge.psd import _current_pulse_model as current_pulse_model

from . import SimflowConfig
from . import metadata as mutils


def _merge_defaults(user: dict, default: dict) -> dict:
    # merge default into user without overwriting user values
    result = dict(default)
    for k, v in user.items():
        if k in result and isinstance(result[k], dict) and isinstance(v, dict):
            result[k] = _merge_defaults(v, result[k])
        else:
            result[k] = v
    return result


def init_simflow_context(raw_config: dict, workflow) -> AttrsDict:
    """Pre-process and sanitize the Simflow configuration.

    - set default configuration fields;
    - substitute ``$_`` and environment variables;
    - convert to :class:`~dbetto.attrsdict.AttrsDict`;
    - cast filesystem paths to :class:`pathlib.Path`;
    - clone and configure `legend-metadata`;
    - attach a :class:`legendmeta.LegendMetadata` instance to the Simflow
      configuration;
    - export important environment variables.

    Returns a dictionary with useful objects to be used in the Simflow
    Snakefiles (i.e. the "context").
    """
    if not raw_config:
        msg = "you must set a config file with --configfile"
        raise RuntimeError(msg)

    raw_config = _merge_defaults(
        raw_config,
        {"benchmark": {"enabled": False}, "nersc": {"dvs_ro": False, "scratch": False}},
    )

    ldfs.workflow.utils.subst_vars_in_snakemake_config(workflow, raw_config)
    config = AttrsDict(raw_config)

    # convert all strings in the "paths" block to pathlib.Path
    def _make_path(d):
        for k, v in d.items():
            if isinstance(v, str):
                d[k] = Path(v)
            else:
                d[k] = _make_path(v)
        return d

    config["paths"] = _make_path(config.paths)

    # NOTE: this will attempt a clone of legend-metadata, if the directory does not exist
    # NOTE: don't use lazy=True, we need a fully functional TextDB
    metadata = LegendMetadata(config.paths.metadata)
    if "legend_metadata_version" in config:
        metadata.checkout(config.legend_metadata_version)
    config["metadata"] = metadata

    # make sure all simflow plots are made with a consistent style
    # I have verified only that this variable is visible in scripts (not shell directives)
    os.environ["MPLCONFIGDIR"] = f"{workflow.basedir}/src/legendsimflow"

    proctime = (
        "benchmark"
        if config.benchmark.enabled
        else datetime.now().strftime("%Y%m%dT%H%M%SZ")
    )

    return AttrsDict(
        {
            "config": config,
            "basedir": workflow.basedir,
            "proctime": proctime,
        }
    )


def _get_matching_volumes(volume_list: list, patterns: str | list) -> list[str]:
    """Get the list of volumes from the GDML. The string can include wildcards."""

    wildcard_list = [patterns] if isinstance(patterns, str) else patterns

    # find all volumes matching at least one pattern
    matched_list = []
    for key in volume_list:
        for name in wildcard_list:
            if fnmatch.fnmatch(key, name):
                matched_list.append(key)
    return matched_list


def get_lar_minishroud_confine_commands(
    reg: pg4.geant4.Registry,
    pattern: str = "minishroud_side*",
    inside: bool = True,
    lar_name: str = "LAr",
) -> list[str]:
    """Extract the commands for the LAr confinement inside/ outside the NMS from the GDML.

    Parameters
    ----------
    reg
        The registry describing the geometry.
    pattern
        The pattern used to search for physical volumes of minishrouds.
    inside
        If True, generate points inside the minishroud (NMS) volumes; if False,
        exclude the minishroud volumes from the generation region.
    lar_name
        The name of the physical volume of the LAr.

    Returns
    -------
    a list of confinement commands for remage.
    """
    string_list = _get_matching_volumes(list(reg.physicalVolumeDict.keys()), pattern)
    # correct sampling mode
    mode = "IntersectPhysicalWithGeometrical" if inside else "SubtractGeometrical"

    # physical volume sampling
    lines = [f"/RMG/Generator/Confinement/SamplingMode {mode}"]
    lines += [f"/RMG/Generator/Confinement/Physical/AddVolume {lar_name}"]

    for s in string_list:
        vol = reg.physicalVolumeDict[s]

        center = vol.position.eval()
        solid = vol.logicalVolume.solid

        # Validate expected geometry structure before accessing attributes
        if not hasattr(solid, "obj1") or solid.obj1 is None:
            raise ValueError(
                f"Expected solid for physical volume '{s}' to have an 'obj1' "
                "attribute representing the outer minishroud cylinder, but it was "
                f\"{'missing' if not hasattr(solid, 'obj1') else 'None'}\".
            )

        outer_ms = solid.obj1

        if not hasattr(outer_ms, "pRMax") or not hasattr(outer_ms, "pDz"):
            missing = [
                name for name in ("pRMax", "pDz") if not hasattr(outer_ms, name)
            ]
            raise ValueError(
                "Outer minishroud solid for physical volume "
                f"'{s}' is missing expected attribute(s): {', '.join(missing)}. "
                "The confinement geometry assumes a cylindrical solid with "
                "'pRMax' and 'pDz' attributes."
            )
        r_max = outer_ms.pRMax
        dz = outer_ms.pDz

        # type conversions from pyg4ometry types
        if not isinstance(r_max, float):
            r_max = r_max.eval()

        if not isinstance(dz, float):
            dz = dz.eval()

        command = "AddSolid" if inside else "AddExcludeSolid"
        lines.append(f"/RMG/Generator/Confinement/Geometrical/{command} Cylinder")

        lines.append(
            f"/RMG/Generator/Confinement/Geometrical/CenterPositionX {center[0]} mm"
        )
        lines.append(
            f"/RMG/Generator/Confinement/Geometrical/CenterPositionY {center[1]} mm"
        )
        lines.append(
            f"/RMG/Generator/Confinement/Geometrical/CenterPositionZ {center[2]} mm"
        )
        lines.append(
            f"/RMG/Generator/Confinement/Geometrical/Cylinder/OuterRadius {r_max} mm"
        )
        lines.append(
            f"/RMG/Generator/Confinement/Geometrical/Cylinder/Height {2 * dz} mm"
        )

    return lines


def setup_logdir_link(config: SimflowConfig, proctime):
    """Set up the timestamp-tagged directory for the workflow log files."""
    logdir = Path(config.paths.log)
    logdir.mkdir(parents=True, exist_ok=True)

    # create a handy link to access latest log directory
    link = logdir / "latest"
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(proctime, target_is_directory=True)


# FIXME: this should be removed once the PRL25 data is reprocessed
def _get_lh5_table(
    metadata: LegendMetadata,
    fname: str | Path,
    hpge: str,
    tier: str,
    runid: str,
) -> str:
    """The correct LH5 table path.

    Determines the correct path to a `hpge` detector table in tier `tier`.
    """
    # check if the latest format is available
    path = f"{tier}/{hpge}"
    if lh5.ls(fname, path) == [path]:
        return path

    # otherwise fall back to the old format
    timestamp = mutils.runinfo(metadata, runid).start_key
    rawid = metadata.channelmap(timestamp)[hpge].daq.rawid
    return f"ch{rawid}/{tier}"


def _curve_fit_popt_to_dict(popt: ArrayLike) -> dict:
    """Get the :func:`scipy.curve_fit` parameter results as a dictionary"""
    params = list(inspect.signature(current_pulse_model).parameters)
    param_names = params[1:]

    popt_dict = dict(zip(param_names, popt, strict=True))
    popt_dict["mean_aoe"] = popt_dict["amax"] / 1593

    for key, value in popt_dict.items():
        popt_dict[key] = float(f"{value:.3g}")

    return popt_dict
