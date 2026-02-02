from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest
from dbetto import AttrsDict
from lgdo import Array, Table

from legendsimflow import utils


def test_get_parameter_dict():
    popt = [100, 10, 60, 0.6, 100, 0.2, 60]

    popt_dict = utils._curve_fit_popt_to_dict(popt)
    assert len(popt_dict) == 7


def test_hash_string_int():
    int_hash = utils.string_to_remage_seed("blah blah legend simflow")
    assert isinstance(int_hash, int)
    assert int_hash >= 0


def test_add_field_string():
    tab = Table(size=10)

    utils.add_field_string("string", tab, "luigi was here")

    assert isinstance(tab.string, Array)
    assert np.all(tab.string.view_as("np").astype("str") == "luigi was here")


def test_get_dataflow_config(test_l200data):
    config = utils.lookup_dataflow_config(test_l200data / "v2.1.5")

    assert isinstance(config, AttrsDict)
    assert "paths" in config
    assert "$_" not in config.paths.tier


def test_check_nans_leq():
    utils.check_nans_leq([1, 2, 3, 4], "boh", 0.1)
    utils.check_nans_leq([[1], [2, 3], 4], "boh", 0.1)
    utils.check_nans_leq([[1], [2, 3], np.nan], "boh", 0.5)

    with pytest.raises(RuntimeError):
        utils.check_nans_leq([[1], [2, 3], np.nan], "boh", 0.1)


def test_merge_defaults():
    """Test _merge_defaults function for recursive dictionary merging."""
    # Test basic merge
    user = {"a": 1, "b": 2}
    default = {"b": 999, "c": 3}
    result = utils._merge_defaults(user, default)
    assert result == {"a": 1, "b": 2, "c": 3}

    # Test nested merge
    user = {"outer": {"inner": 1, "keep": "user"}}
    default = {"outer": {"inner": 999, "default_key": "default"}, "top": "level"}
    result = utils._merge_defaults(user, default)
    assert result["outer"]["inner"] == 1
    assert result["outer"]["keep"] == "user"
    assert result["outer"]["default_key"] == "default"
    assert result["top"] == "level"

    # Test empty dictionaries
    assert utils._merge_defaults({}, {"a": 1}) == {"a": 1}
    assert utils._merge_defaults({"a": 1}, {}) == {"a": 1}

    # Test overwriting non-dict with dict
    user = {"key": {"nested": "value"}}
    default = {"key": "simple"}
    result = utils._merge_defaults(user, default)
    assert result["key"] == {"nested": "value"}


def test_setup_logdir_link():
    """Test setup_logdir_link creates symlink to log directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        log_path = tmpdir_path / "logs"

        # Create a mock config object
        config = AttrsDict({"paths": {"log": log_path}})
        proctime = "20240101T120000Z"

        # Create the target directory
        log_path.mkdir(parents=True, exist_ok=True)
        (log_path / proctime).mkdir(parents=True, exist_ok=True)

        # First call should create the symlink
        utils.setup_logdir_link(config, proctime)

        link = log_path / "latest"
        assert link.exists()
        assert link.is_symlink()
        assert link.readlink() == Path(proctime)

        # Second call should update the symlink
        new_proctime = "20240102T130000Z"
        (log_path / new_proctime).mkdir(parents=True, exist_ok=True)
        utils.setup_logdir_link(config, new_proctime)

        assert link.exists()
        assert link.is_symlink()
        assert link.readlink() == Path(new_proctime)


def test_hash_dict():
    """Test hash_dict produces consistent JSON representation."""
    # Test basic dictionary
    d1 = {"a": 1, "b": 2, "c": 3}
    hash1 = utils.hash_dict(d1)
    assert isinstance(hash1, str)
    assert '"a": 1' in hash1
    assert '"b": 2' in hash1

    # Test order independence
    d2 = {"c": 3, "a": 1, "b": 2}
    hash2 = utils.hash_dict(d2)
    assert hash1 == hash2

    # Test nested dictionary
    d3 = {"outer": {"inner": 1}, "list": [1, 2, 3]}
    hash3 = utils.hash_dict(d3)
    assert isinstance(hash3, str)
    assert "outer" in hash3

    # Test AttrsDict
    attrs_dict = AttrsDict({"key": "value", "num": 42})
    hash4 = utils.hash_dict(attrs_dict)
    assert isinstance(hash4, str)
    assert "key" in hash4


def test_get_hit_tier_name(test_l200data):
    """Test get_hit_tier_name returns correct tier."""
    # The test data config has tier_hit and tier_pht defined
    # but the actual directories don't exist, so we expect a RuntimeError
    with pytest.raises(RuntimeError, match="does not contain a valid pht or hit tier"):
        utils.get_hit_tier_name(str(test_l200data / "v2.1.5"))
