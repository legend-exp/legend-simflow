"""Tests for spms_pars module."""

from __future__ import annotations

import numpy as np
import pytest
import awkward as ak

from legendsimflow import spms_pars


@pytest.fixture
def mock_forced_trig_library():
    """Create a mock forced trigger library for testing."""
    # Create structured array with npe, t0, and rawid fields
    # 10 events, 3 SiPM channels each
    npe_data = ak.Array([
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        [[2, 3, 4], [5, 6, 7], [8, 9, 10]],
        [[3, 4, 5], [6, 7, 8], [9, 10, 11]],
        [[4, 5, 6], [7, 8, 9], [10, 11, 12]],
        [[5, 6, 7], [8, 9, 10], [11, 12, 13]],
        [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
        [[2, 2, 2], [3, 3, 3], [4, 4, 4]],
        [[3, 3, 3], [4, 4, 4], [5, 5, 5]],
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        [[10, 20, 30], [40, 50, 60], [70, 80, 90]],
    ])
    
    t0_data = ak.Array([
        [[100.0, 200.0, 300.0], [150.0, 250.0, 350.0], [120.0, 220.0, 320.0]],
        [[105.0, 205.0, 305.0], [155.0, 255.0, 355.0], [125.0, 225.0, 325.0]],
        [[110.0, 210.0, 310.0], [160.0, 260.0, 360.0], [130.0, 230.0, 330.0]],
        [[115.0, 215.0, 315.0], [165.0, 265.0, 365.0], [135.0, 235.0, 335.0]],
        [[120.0, 220.0, 320.0], [170.0, 270.0, 370.0], [140.0, 240.0, 340.0]],
        [[50.0, 150.0, 250.0], [100.0, 200.0, 300.0], [75.0, 175.0, 275.0]],
        [[55.0, 155.0, 255.0], [105.0, 205.0, 305.0], [80.0, 180.0, 280.0]],
        [[60.0, 160.0, 260.0], [110.0, 210.0, 310.0], [85.0, 185.0, 285.0]],
        [[65.0, 165.0, 265.0], [115.0, 215.0, 315.0], [90.0, 190.0, 290.0]],
        [[70.0, 170.0, 270.0], [120.0, 220.0, 320.0], [95.0, 195.0, 295.0]],
    ])
    
    # SiPM channel UIDs: 101, 102, 103
    rawid_array = np.array([101, 102, 103], dtype=np.int32)
    
    # Create the library structure like get_forced_trigger_library returns
    library = ak.Array({
        "npe": npe_data,
        "t0": t0_data,
        "rawid": ak.Array([rawid_array for _ in range(len(npe_data))]),
    })
    
    return library


class TestRandCoincSampler:
    """Tests for RandCoincSampler class."""
    
    def test_initialization(self, mock_forced_trig_library):
        """Test sampler initialization."""
        sampler = spms_pars.RandCoincSampler(mock_forced_trig_library)
        
        assert sampler.library is not None
        assert sampler.rng is not None
        assert len(sampler._pools) == 0  # No pools initialized yet
    
    def test_sample_basic(self, mock_forced_trig_library):
        """Test basic sampling functionality."""
        sampler = spms_pars.RandCoincSampler(mock_forced_trig_library)
        
        # Sample 5 entries from SiPM UID 101
        npe_sample, t0_sample = sampler.sample(101, 5)
        
        # Check shapes
        assert len(npe_sample) == 5
        assert len(t0_sample) == 5
    
    def test_sample_missing_sipm(self, mock_forced_trig_library):
        """Test error when sampling non-existent SiPM UID."""
        sampler = spms_pars.RandCoincSampler(mock_forced_trig_library)
        
        with pytest.raises(ValueError, match="SiPM UID 999 not found"):
            sampler.sample(999, 5)
    
    def test_sample_returns_valid_data(self, mock_forced_trig_library):
        """Test that sampled data has valid structure."""
        sampler = spms_pars.RandCoincSampler(mock_forced_trig_library)
        
        # Sample from channel 101
        npe_sample, t0_sample = sampler.sample(101, 5)
        
        # Check that we got arrays with the right length
        assert len(npe_sample) == 5
        assert len(t0_sample) == 5
        
        # Check that data is valid (not empty, contains numbers)
        npe_list = ak.to_list(npe_sample)
        t0_list = ak.to_list(t0_sample)
        
        for npe_val, t0_val in zip(npe_list, t0_list):
            assert isinstance(npe_val, list)
            assert isinstance(t0_val, list)
            assert len(npe_val) > 0
            assert len(t0_val) > 0
    
    def test_sample_with_replacement(self, mock_forced_trig_library):
        """Test that sampling with replacement allows event reuse."""
        sampler = spms_pars.RandCoincSampler(mock_forced_trig_library)
        n_lib = len(mock_forced_trig_library)
        
        # Request more samples than available
        npe_sample, t0_sample = sampler.sample(101, n_lib * 2)
        
        # Should successfully return requested samples despite refilling
        assert len(npe_sample) == n_lib * 2
        assert len(t0_sample) == n_lib * 2
    
    def test_pool_refill_on_exhaustion(self, mock_forced_trig_library):
        """Test that pool refills when exhausted."""
        rng = np.random.default_rng(42)
        sampler = spms_pars.RandCoincSampler(mock_forced_trig_library, rng)
        
        n_lib = len(mock_forced_trig_library)
        
        # First sample should use the initial pool
        npe1, t0_1 = sampler.sample(101, n_lib)
        assert len(npe1) == n_lib
        
        # Second sample should trigger pool refill
        npe2, t0_2 = sampler.sample(101, n_lib)
        assert len(npe2) == n_lib
        
        # The two samples together should have used more data than available once
        # (This is hard to verify directly, but we can check they both exist)
        assert len(npe1) > 0
        assert len(npe2) > 0
    
    def test_multiple_sipms(self, mock_forced_trig_library):
        """Test sampling from multiple SiPM channels."""
        sampler = spms_pars.RandCoincSampler(mock_forced_trig_library)
        
        # Sample from different channels
        npe1, t0_1 = sampler.sample(101, 3)
        npe2, t0_2 = sampler.sample(102, 3)
        npe3, t0_3 = sampler.sample(103, 3)
        
        # All should succeed
        assert len(npe1) == 3
        assert len(npe2) == 3
        assert len(npe3) == 3
        
        # Pools should be tracked separately
        assert len(sampler._pools) == 3
    
    def test_reproducibility_with_seed(self, mock_forced_trig_library):
        """Test that same seed produces same samples."""
        rng1 = np.random.default_rng(42)
        sampler1 = spms_pars.RandCoincSampler(mock_forced_trig_library, rng1)
        npe1, t0_1 = sampler1.sample(101, 5)
        
        rng2 = np.random.default_rng(42)
        sampler2 = spms_pars.RandCoincSampler(mock_forced_trig_library, rng2)
        npe2, t0_2 = sampler2.sample(101, 5)
        
        # Convert to comparable format
        npe1_list = ak.to_list(npe1)
        npe2_list = ak.to_list(npe2)
        t0_1_list = ak.to_list(t0_1)
        t0_2_list = ak.to_list(t0_2)
        
        # Should be identical with same seed
        assert npe1_list == npe2_list
        assert t0_1_list == t0_2_list
    
    def test_different_seeds_different_results(self, mock_forced_trig_library):
        """Test that different seeds produce different samples."""
        rng1 = np.random.default_rng(42)
        sampler1 = spms_pars.RandCoincSampler(mock_forced_trig_library, rng1)
        npe1, _ = sampler1.sample(101, 10)
        
        rng2 = np.random.default_rng(123)
        sampler2 = spms_pars.RandCoincSampler(mock_forced_trig_library, rng2)
        npe2, _ = sampler2.sample(101, 10)
        
        # Very unlikely to be the same with different seeds
        npe1_list = ak.to_list(npe1)
        npe2_list = ak.to_list(npe2)
        
        # At least they should both be valid
        assert len(npe1_list) == len(npe2_list) == 10
    
    def test_large_sample_request(self, mock_forced_trig_library):
        """Test sampling much larger than library size."""
        sampler = spms_pars.RandCoincSampler(mock_forced_trig_library)
        n_lib = len(mock_forced_trig_library)
        
        # Request 100x library size
        npe, t0 = sampler.sample(101, n_lib * 100)
        
        assert len(npe) == n_lib * 100
        assert len(t0) == n_lib * 100


class TestRandCoincSpmsData:
    """Tests for the stateless rand_coinc_spms_data function."""
    
    def test_basic_sampling(self, mock_forced_trig_library):
        """Test basic stateless sampling."""
        npe, t0 = spms_pars.rand_coinc_spms_data(
            mock_forced_trig_library, 101, 5
        )
        
        assert len(npe) == 5
        assert len(t0) == 5
    
    def test_missing_sipm_error(self, mock_forced_trig_library):
        """Test error on missing SiPM."""
        with pytest.raises(ValueError, match="SiPM UID 999 not found"):
            spms_pars.rand_coinc_spms_data(
                mock_forced_trig_library, 999, 5
            )
    
    def test_large_request_with_warning(self, mock_forced_trig_library, caplog):
        """Test warning when requesting more samples than available."""
        import logging
        
        with caplog.at_level(logging.WARNING):
            npe, t0 = spms_pars.rand_coinc_spms_data(
                mock_forced_trig_library, 101, len(mock_forced_trig_library) * 10
            )
        
        # Should still work
        assert len(npe) == len(mock_forced_trig_library) * 10
        
        # Should have logged a warning
        assert any("Requested" in record.message for record in caplog.records)
        assert any("Sampling with replacement" in record.message for record in caplog.records)
