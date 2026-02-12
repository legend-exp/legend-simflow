from __future__ import annotations

from pathlib import Path
import re
import logging
import awkward as ak
import numpy as np

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
    
    evt_path = Path(
        df_cfg[f"tier_{evt_tier_name}"]
    ).resolve()
    evt_files = list((evt_path / data_type / period / run).glob("*"))
    
    return evt_files


class RandCoincSampler:
    """Stateful sampler for random coincidence SiPM data.
    
    Minimizes reuse of events across multiple sampling calls by maintaining
    a pool of unused indices for each SiPM channel.
    """
    
    def __init__(self, forced_trig_library: ak.Array, rng: np.random.Generator = None):
        """Initialize the sampler.
        
        Parameters
        ----------
        forced_trig_library
            Library of forced trigger events containing npe, t0, and rawid fields.
        rng
            Random number generator. If None, uses default_rng().
        """
        self.library = forced_trig_library
        self.rng = rng if rng is not None else np.random.default_rng()
        
        # Map from sipm_uid to (channel_index, available_indices)
        self._pools = {}
        
    def _get_pool(self, sipm_uid: int) -> tuple[int, list[int]]:
        """Get or create the index pool for a SiPM channel."""
        if sipm_uid not in self._pools:
            # Find the channel index for this SiPM UID
            channel_indices = ak.where(self.library.rawid[0] == sipm_uid)[0]
            
            if len(channel_indices) == 0:
                msg = f"SiPM UID {sipm_uid} not found in forced trigger library"
                raise ValueError(msg)
            
            ch_idx = int(channel_indices[0])
            
            # Initialize with shuffled indices
            n_available = len(self.library.npe[:, ch_idx])
            if n_available == 0:
                msg = f"No data available for SiPM UID {sipm_uid}"
                raise ValueError(msg)
            
            indices = list(range(n_available))
            self.rng.shuffle(indices)
            
            self._pools[sipm_uid] = (ch_idx, indices)
        
        return self._pools[sipm_uid]
    
    def sample(self, sipm_uid: int, n_entries: int) -> tuple[ak.Array, ak.Array]:
        """Sample random coincidence data for a SiPM channel.
        
        Draws from unused events first. When the pool is exhausted, it refills
        with all available events (shuffled) and continues sampling.
        
        Parameters
        ----------
        sipm_uid
            SiPM channel ID to sample.
        n_entries
            Number of entries to sample.
        
        Returns
        -------
        npe_sample
            Sampled photoelectron counts.
        t0_sample
            Sampled times.
        """
        ch_idx, pool = self._get_pool(sipm_uid)
        
        # Select data for this channel
        npe = self.library.npe[:, ch_idx]
        t0 = self.library.t0[:, ch_idx]
        n_available = len(npe)
        
        sampled_indices = []
        
        while len(sampled_indices) < n_entries:
            # If pool is empty, refill it
            if len(pool) == 0:
                log.debug(
                    f"Refilling pool for SiPM UID {sipm_uid} "
                    f"(need {n_entries - len(sampled_indices)} more samples)"
                )
                pool[:] = list(range(n_available))
                self.rng.shuffle(pool)
            
            # Draw as many as we need (or as many as available)
            n_to_draw = min(len(pool), n_entries - len(sampled_indices))
            sampled_indices.extend(pool[:n_to_draw])
            del pool[:n_to_draw]
        
        # Convert to array and sample
        indices = np.array(sampled_indices)
        npe_sample = npe[indices]
        t0_sample = t0[indices]
        
        return npe_sample, t0_sample


def rand_coinc_spms_data(
    forced_trig_library: ak.Array,
    sipm_uid: int,
    n_entries: int,
    rng: np.random.Generator = None,
) -> tuple[ak.Array, ak.Array]:
    """Sample random coincidence SiPM data from forced trigger library.
    
    Note: This is a stateless convenience function. For processing data in chunks
    where you want to minimize event reuse across chunks, use RandCoincSampler instead.
    
    Parameters
    ----------
    forced_trig_library
        Library of forced trigger events containing npe, t0, and rawid fields.
    sipm_uid
        SiPM channel ID to sample.
    n_entries
        Number of entries to sample.
    rng
        Random number generator. If None, uses default_rng().
    
    Returns
    -------
    npe_sample
        Sampled photoelectron counts.
    t0_sample
        Sampled times.
    """
    if rng is None:
        rng = np.random.default_rng()
    
    # Find the channel index for this SiPM UID
    # rawid[0] gives the channel IDs (should be the same for all events)
    channel_indices = ak.where(forced_trig_library.rawid[0] == sipm_uid)[0]
    
    if len(channel_indices) == 0:
        msg = f"SiPM UID {sipm_uid} not found in forced trigger library"
        raise ValueError(msg)
    
    ch_idx = int(channel_indices[0])
    
    # Select data for this channel from all events
    npe = forced_trig_library.npe[:, ch_idx]
    t0 = forced_trig_library.t0[:, ch_idx]
    
    # Sample with replacement
    n_available = len(npe)
    if n_available == 0:
        msg = f"No data available for SiPM UID {sipm_uid}"
        raise ValueError(msg)
    
    if n_entries > n_available:
        log.warning(
            f"Requested {n_entries} samples but only {n_available} available for SiPM UID {sipm_uid}. "
            "Sampling with replacement - some events will be reused."
        )
    
    indices = rng.choice(n_available, size=n_entries, replace=True)
    
    npe_sample = npe[indices]
    t0_sample = t0[indices]
    
    return npe_sample, t0_sample
    
    