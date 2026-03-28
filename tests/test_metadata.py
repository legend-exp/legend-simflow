from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
from dbetto import AttrsDict

from legendsimflow import metadata
from legendsimflow.exceptions import SimflowConfigError


def test_all(config):
    assert isinstance(
        metadata.get_simconfig(config, "stp", "birds_nest_K40"), AttrsDict
    )
    assert (
        metadata.get_simconfig(config, "stp", "birds_nest_K40", "number_of_jobs") == 2
    )
    assert isinstance(
        metadata.smk_hash_simconfig(config, None, tier="stp", simid="birds_nest_K40"),
        str,
    )

    assert "livetime_in_s" in metadata.runinfo(config.metadata, "l200-p02-r000-phy")

    assert (
        "operational_voltage_in_V"
        in metadata.simpars(config.metadata, "geds.opv", "l200-p02-r002-phy").V99000A
    )

    assert isinstance(
        metadata.get_vtx_simconfig(config, "lar_hpge_shell_K42"), AttrsDict
    )


def test_run_stuff(config):
    assert metadata.parse_runid("l200-p42-r999-ant") == ("l200", 42, 999, "ant")

    assert metadata.is_runid("l200-p00-r000-phy")
    assert not metadata.is_runid("l200-p00-r00-phy")
    assert not metadata.is_runid("l200-p00-r000-ph0")
    assert not metadata.is_runid("l200-k00-r000-phy")
    assert not metadata.is_runid("l1000-p00-r000-phy")

    assert metadata.query_runlist_db(config.metadata, "valid.phy.p02") == [
        f"l200-p02-r00{r}-phy" for r in range(8)
    ]
    assert metadata.expand_runlist(config.metadata, "~runlists:valid.phy.p02") == [
        f"l200-p02-r00{r}-phy" for r in range(8)
    ]
    assert metadata.expand_runlist(
        config.metadata, ["~runlists:valid.phy.p02", "~runlists:valid.cal.p02"]
    ) == sorted([f"l200-p02-r00{r}-{dt}" for r in range(8) for dt in ("phy", "cal")])

    assert metadata.get_runlist(config, "exotic_physics_hpge") == config.runlist
    assert metadata.get_runlist(config, "phbr_surface_Ra228_to_Ac228") == [
        "l200-p02-r000-phy",
        "l200-p02-r001-phy",
    ]

    assert metadata.get_runlist(config, "sis1_z8224_src1_Ra224_to_Pb208") == [
        f"l200-p02-r00{r}-phy" for r in range(8)
    ]

    assert (
        metadata.reference_cal_run(config.metadata, "l200-p16-r006-phy")
        == "l200-p16-r006-cal"
    )
    assert (
        metadata.reference_cal_run(config.metadata, "l200-p16-r008-ssc")
        == "l200-p16-r007-cal"
    )
    assert (
        metadata.reference_cal_run(config.metadata, "l200-p16-r009-ssc")
        == "l200-p16-r007-cal"
    )


def test_is_simid():
    # valid simids
    assert metadata.is_simid("hpge_bulk_Rn222_to_Po214")
    assert metadata.is_simid("birds_nest_K40")
    assert metadata.is_simid("simid123")
    assert metadata.is_simid("with-hyphen")
    assert metadata.is_simid("a")

    # invalid simids
    assert not metadata.is_simid("has.dot")
    assert not metadata.is_simid("has dot")
    assert not metadata.is_simid("has@special")
    assert not metadata.is_simid("")
    assert not metadata.is_simid("tier.simid")


def test_validate_simconfig_keys():
    # all valid keys — should not raise
    valid = {
        "hpge_bulk_Rn222_to_Po214": {},
        "birds_nest_K40": {},
        "simid-with-hyphens": {},
    }
    metadata.validate_simconfig_keys(valid)

    # one invalid key (contains a dot)
    invalid = {"valid_key": {}, "bad.key": {}}
    with pytest.raises(SimflowConfigError, match=r"bad\.key"):
        metadata.validate_simconfig_keys(invalid)

    # multiple invalid keys, block label included in message
    multi_invalid = {"ok": {}, "also.bad": {}, "has space": {}}
    with pytest.raises(SimflowConfigError, match=r"also\.bad"):
        metadata.validate_simconfig_keys(multi_invalid, block="test.block")


def test_get_simconfig_validates_keys(config):
    # getting the full stp simconfig should succeed (all keys are valid simids)
    simcfg = metadata.get_simconfig(config, "stp")
    assert all(metadata.is_simid(k) for k in simcfg)


def test_get_simconfig_invalid_key_raises(config):
    # inject an invalid top-level key and ensure full simconfig load validates it
    config.metadata.simprod.config.tier.stp[config.experiment].simconfig[123] = {}

    with pytest.raises(SimflowConfigError, match=r"123"):
        metadata.get_simconfig(config, "stp")


def test_encode_usability():
    assert metadata.encode_usability("on") == 0
    assert metadata.encode_usability("ac") == 1
    assert metadata.encode_usability("off") == 2

    for use in ["on", "off", "ac"]:
        assert metadata.decode_usability(metadata.encode_usability(use)) == use


def test_fccd(config):
    assert metadata.get_sanitized_fccd(config.metadata, "B99000A") == 0.75


def test_extract_integer():
    """Test extract_integer reads integer from file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        test_file = tmpdir_path / "test_int.txt"

        # Test simple integer
        test_file.write_text("42")
        assert metadata.extract_integer(test_file) == 42

        # Test with whitespace
        test_file.write_text("  123  \n")
        assert metadata.extract_integer(test_file) == 123

        # Test negative integer
        test_file.write_text("-999")
        assert metadata.extract_integer(test_file) == -999
