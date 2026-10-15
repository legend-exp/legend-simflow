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

"""Checks of the Simflow configuration against the LEGEND-200 data production.

The data production is found at ``paths.l200data`` in the Simflow
configuration. Each check does nothing if this path is not set.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path

from dbetto import AttrsDict
from legendmeta import LegendMetadata

from . import SimflowConfig, aggregate
from .exceptions import SimflowConfigError

log = logging.getLogger(__name__)


def detinfo_against_l200data(
    config: SimflowConfig,
    detinfo: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> None:
    """Check that the detector usabilities match the data production metadata.

    Does nothing if ``paths.l200data`` is not set. Logs a warning instead of
    raising an error if ``allow_metadata_differences`` is true.

    Parameters
    ----------
    config
        Simflow configuration object.
    detinfo
        Detector info computed from the Simflow metadata, as returned by
        :func:`.aggregate.pivot_detinfo`.
    """
    if "l200data" not in config.paths:
        return

    inputs = Path(config.paths.l200data) / "inputs"
    data_config = dict(config)
    data_config["metadata"] = LegendMetadata(str(inputs), lazy=True)
    data_config = AttrsDict(data_config)

    flags = ["usability", "psd_usability"]
    # the data production metadata has no simprod config: query the same runs
    runids = {runid for flag in flags for runid in detinfo.get(flag, {})}
    data_detinfo = aggregate.pivot_detinfo(
        aggregate.gen_list_of_all_usabilities(data_config, runids).to_dict()
    )

    for flag in flags:
        sim = detinfo.get(flag, {})
        data = data_detinfo.get(flag, {})
        if sim == data:
            continue

        diffs = []
        for runid in sorted(set(sim) | set(data)):
            sim_run = sim.get(runid, {})
            data_run = data.get(runid, {})
            for det in sorted(set(sim_run) | set(data_run)):
                sim_val = sim_run.get(det, "<missing>")
                data_val = data_run.get(det, "<missing>")
                if sim_val != data_val:
                    diffs.append(
                        f"  {runid} {det}: Simflow metadata has {sim_val!r}, "
                        f"data production has {data_val!r}"
                    )

        msg = (
            f"'{flag}' differs between the Simflow metadata and the data "
            f"production metadata in {inputs} for {len(diffs)} detector(s):\n"
            + "\n".join(diffs)
        )
        if config.get("allow_metadata_differences", False):
            log.warning(msg)
        else:
            raise SimflowConfigError(msg)
