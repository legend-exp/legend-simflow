from __future__ import annotations

from pathlib import Path

import yaml

from legendsimflow import utils
from legendsimflow.warmup import _resolve_dsp_config, warm_hpge_dsp_cache

l200data = Path(__file__).parent / "l200data" / "v3.0.0"


def test_lookup_dsp_config():
    cfg = utils.lookup_dsp_config(l200data)
    assert cfg.exists()
    assert "ICPC-dsp_proc_chain" in cfg.name


def test_warm_hpge_dsp_cache():
    """Warm the real current-pulse DSP chain on a small synthetic waveform (must not raise)."""
    warm_hpge_dsp_cache(utils.lookup_dsp_config(l200data))


def test_resolve_dsp_config_skips_without_l200data(tmp_path):
    """No DSP warmup when the config has no ``l200data``."""
    cfg = tmp_path / "simflow-config.yaml"
    cfg.write_text(yaml.safe_dump({"experiment": "l200cfg09", "paths": {}}))
    assert _resolve_dsp_config(cfg) is None


def test_resolve_dsp_config_skips_when_missing():
    """No DSP warmup when the config path is missing or ``None``."""
    assert _resolve_dsp_config("does-not-exist.yaml") is None
    assert _resolve_dsp_config(None) is None
