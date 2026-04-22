from __future__ import annotations

import shutil

import pytest
from lgdo import lh5


@pytest.mark.needs_julia
@pytest.mark.skipif(shutil.which("julia") is None, reason="julia not installed")
def test_make_hpge_drift_time_maps_l1000(legend_dtmap_path):
    assert legend_dtmap_path.exists(), "Drift time map LH5 file was not created"

    top_keys = lh5.ls(legend_dtmap_path)
    assert "V05261B" in top_keys, f"Expected group 'V05261B' in LH5, got: {top_keys}"

    inner_keys = lh5.ls(legend_dtmap_path, "V05261B/")
    for expected in ("drift_time_000_deg", "drift_time_045_deg", "r", "z"):
        assert f"V05261B/{expected}" in inner_keys, (
            f"Expected dataset 'V05261B/{expected}' in LH5, found: {inner_keys}"
        )
