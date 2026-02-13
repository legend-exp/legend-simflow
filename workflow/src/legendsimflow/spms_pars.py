from __future__ import annotations

import logging
import re
from pathlib import Path

import awkward as ak

from . import utils

log = logging.getLogger(__name__)


def lookup_evt_files(l200data: str, runid: str, evt_tier_name: str) -> list[str | Path]:
    """Lookup the paths to the `evt` files."""

    _, period, run, data_type = re.split(r"\W+", runid)

    if isinstance(l200data, str):
        l200data = Path(l200data)

    dataflow_config = utils.lookup_dataflow_config(l200data)

    # get the paths to hit and raw tier files
    df_cfg = (
        dataflow_config["setups"]["l200"]["paths"]
        if ("setups" in dataflow_config)
        else dataflow_config["paths"]
    )

    evt_path = Path(df_cfg[f"tier_{evt_tier_name}"]).resolve()
    return list((evt_path / data_type / period / run).glob("*"))


def forced_trig_sipm_data(
    ft_library: ak.Array,
    sipm: str,
    sipm_uid: int,
) -> tuple[ak.Array, ak.Array]:
    """Extract npe and times for a specific SiPM channel from the library.

    Filters the forced trigger library to return only the photoelectron counts
    and times for a single SiPM channel across all events.

    Parameters
    ----------
    ft_library
        Library of forced trigger events containing npe, t0, and rawid fields.
    sipm
        SiPM channel name to extract data for. If "all", flatten channel dimension.
    sipm_uid
        SiPM channel ID to extract data for.

    Returns
    -------
    npe
        Photoelectron counts for the SiPM channel across all events.
    t0
        Photoelectron times for the SiPM channel across all events.

    Raises
    ------
    ValueError
        If the SiPM UID is not found in the library.
    """
    if sipm == "all":
        npe = ak.flatten(ft_library.npe, axis=-1)
        t0 = ak.flatten(ft_library.t0, axis=-1)

    else:
        # Find the channel index for this SiPM UID
        # rawid[0] gives the channel IDs (should be the same for all events)
        channel_indices = ak.where(ft_library.rawid[0] == sipm_uid)[0]

        if len(channel_indices) == 0:
            msg = f"SiPM UID {sipm_uid} not found in forced trigger library"
            raise ValueError(msg)

        ch_idx = int(channel_indices[0])

        # Select data for this channel from all events
        npe = ft_library.npe[:, ch_idx]
        t0 = ft_library.t0[:, ch_idx]

    return npe, t0
