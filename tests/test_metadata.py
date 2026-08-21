from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
from dbetto import AttrsDict
from legendmeta import LegendMetadata

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

    assert "livetime_in_s" in metadata.get_runinfo(config, "l200-p02-r000-phy")

    assert (
        "operational_voltage_in_V"
        in metadata.simpars(
            config, "geds.opv", "l200-p02-r002-phy", config.experiment
        ).V99000A
    )

    # default= returns the default when the par directory does not exist
    assert (
        metadata.simpars(
            config,
            "geds.nonexistent",
            "l200-p02-r002-phy",
            config.experiment,
            default=None,
        )
        is None
    )
    assert (
        metadata.simpars(
            config,
            "geds.nonexistent",
            "l200-p02-r002-phy",
            config.experiment,
            default=42,
        )
        == 42
    )

    # without default=, a missing par raises
    with pytest.raises((KeyError, LookupError, FileNotFoundError)):
        metadata.simpars(
            config, "geds.nonexistent", "l200-p02-r002-phy", config.experiment
        )

    assert isinstance(
        metadata.get_vtx_simconfig(config, "lar_hpge_shell_K42"), AttrsDict
    )


def test_run_stuff(config):
    assert metadata.parse_runid("l200-p42-r999-ant") == ("l200", 42, 999, "ant")

    assert metadata.is_runid("l200-p00-r000-phy")
    assert metadata.is_runid("l1000-p00-r000-phy")
    assert not metadata.is_runid("l200-p00-r00-phy")
    assert not metadata.is_runid("l200-p00-r000-ph0")
    assert not metadata.is_runid("l200-k00-r000-phy")
    assert not metadata.is_runid("l200.p00.r000.phy")

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
        metadata.reference_cal_run(config, "l200-p16-r006-phy") == "l200-p16-r006-cal"
    )
    assert (
        metadata.reference_cal_run(config, "l200-p16-r008-ssc") == "l200-p16-r006-cal"
    )
    assert (
        metadata.reference_cal_run(config, "l200-p16-r009-ssc") == "l200-p16-r006-cal"
    )


def test_get_simconfig_empty_file():
    """An empty simconfig.yaml is valid: YAML loads it as None, not as {}.

    The hit tier of an experiment that takes its runlist from ``config.runlist``
    has nothing to put in the file.
    """
    config = AttrsDict(
        {
            "experiment": "l1000dsg01",
            "runlist": ["l1000-p01-r000-phy"],
            "metadata": {
                "simprod": {
                    "config": {"tier": {"hit": {"l1000dsg01": {"simconfig": None}}}}
                }
            },
        }
    )

    assert metadata.get_simconfig(config, "hit") == {}

    with pytest.raises(SimflowConfigError, match="not found"):
        metadata.get_simconfig(config, "hit", "some_simid")

    assert metadata.get_runlist(config, "some_simid") == config.runlist


def test_experiment_prefix():
    assert metadata.experiment_prefix("l200cfg09") == "l200"
    assert metadata.experiment_prefix("l1000dsg01") == "l1000"
    assert metadata.experiment_prefix("legend") == metadata.DEFAULT_RUNID_PREFIX


def test_query_runlist_db_prefix(config):
    assert metadata.query_runlist_db(config.metadata, "valid.phy.p03", "l1000") == [
        f"l1000-p03-r00{r}-phy" for r in range(6)
    ]


def test_get_crystal_name(config):
    diodes = config.metadata.hardware.detectors.germanium.diodes
    assert metadata.get_crystal_name(diodes.V05261B) == "V05261"
    assert metadata.get_crystal_name(diodes.B00000A) == "B99000"

    assert (
        metadata.get_crystal_name(diodes.V05261B)
        in config.metadata.hardware.detectors.germanium.crystals
    )


def test_metadata_overlay(l1000_config, config):
    """The experiment metadata lives in ``paths.config``, and not in the clone."""
    overlay = l1000_config.metadata_overlay
    assert overlay is not None
    assert (
        Path(overlay.__path__)
        == metadata.metadata_overlay_dirname(l1000_config).resolve()
    )

    # legend-metadata describes no LEGEND-1000 experiment
    assert config.get("metadata_overlay") is None


def test_lookup_falls_back_to_the_overlay(l1000_config):
    """The runs, the run lists and the channel map come from the overlay."""
    assert "p99" not in l1000_config.metadata.datasets.runinfo

    rinfo = metadata.get_runinfo(l1000_config, "l1000-p99-r000-phy")
    assert rinfo.start_key == "20000102T000000Z"
    assert rinfo.livetime_in_s > 0

    assert metadata.get_runlist(l1000_config, "some_simid") == [
        "l1000-p99-r000-phy",
        "l1000-p99-r001-phy",
    ]

    chmap = metadata.get_channelmap(l1000_config, rinfo.start_key)
    assert set(chmap.group("system")) == {"geds", "spms"}
    assert chmap.V99900A.analysis.usability == "on"
    # the diode file merged into the channel map comes from the overlay too
    assert chmap.V99900A.production.crystal == "900"


def test_lookup_prefers_legend_metadata(l1000_config):
    """A detector in both databases resolves to the legend-metadata entry."""
    diodes = l1000_config.metadata.hardware.detectors.germanium.diodes
    assert "V05261B" in diodes
    assert metadata.get_diode(l1000_config, "V05261B") == diodes.V05261B

    # V99900A exists only in the overlay
    assert "V99900A" not in diodes
    assert metadata.get_diode(l1000_config, "V99900A").production.crystal == "900"
    assert metadata.get_crystal(l1000_config, "V99900").slices.A.status == "valid"

    with pytest.raises(FileNotFoundError):
        metadata.get_diode(l1000_config, "V00000Z")

    assert metadata.get_diode(l1000_config, "V00000Z", default=None) is None


def test_validate_metadata_overlay_rejects_run_collision(tmp_path, config):
    """Refuse an overlay period that legend-metadata also defines.

    The Simflow looks a run up by period and run number alone. Such a period
    therefore resolves to the legend-metadata run instead.
    """
    overlay_dir = tmp_path / "metadata/l1000dsg01"
    (overlay_dir / "datasets").mkdir(parents=True)
    (overlay_dir / "datasets/runinfo.yaml").write_text(
        "p02:\n  r000:\n    phy:\n      start_key: 20000102T000000Z\n"
    )
    overlay = LegendMetadata(overlay_dir)

    with pytest.raises(SimflowConfigError, match="p02"):
        metadata.validate_metadata_overlay(config.metadata, overlay, overlay_dir)


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
    assert metadata.get_sanitized_fccd(config, "B99000A") == 0.75


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


def test_get_tier_settings_hit(config):
    """get_tier_settings returns the settings object for the hit tier."""
    settings = metadata.get_tier_settings(config, "hit")

    assert "dead_layer_fraction" in settings
    assert "buffer_len" in settings

    assert settings.dead_layer_fraction == 0.5
    assert settings.buffer_len == "500*MB"


def test_get_tier_settings_evt(config):
    """get_tier_settings returns the settings object for the evt tier."""
    settings = metadata.get_tier_settings(config, "evt")

    assert "geds_energy_thr_kev" in settings
    assert "spms_energy_thr_pe" in settings
    assert "buffer_len" in settings

    assert settings.geds_energy_thr_kev == 25
    assert settings.spms_energy_thr_pe == 0
    assert settings.buffer_len == "50*MB"


def test_get_par_settings(config):
    """get_par_settings returns the settings object for the par directory."""
    settings = metadata.get_par_settings(config, "ssd")

    assert isinstance(settings, AttrsDict)
    assert settings.grid_size_in_mm == 10.0
    assert settings.ssd_refinement_limits == [0.2, 0.1, 0.05, 0.02]
    assert settings.padding == 3

    settings_empty = metadata.get_par_settings(config, "nonexistent_par")
    assert isinstance(settings_empty, AttrsDict)
    assert len(settings_empty) == 0
