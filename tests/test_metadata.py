from __future__ import annotations

from dbetto import AttrsDict

from legendsimflow import metadata


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


def test_encode_usability():
    assert metadata.encode_usability("on") == 0
    assert metadata.encode_usability("ac") == 1
    assert metadata.encode_usability("off") == 2

    for use in ["on", "off", "ac"]:
        assert metadata.decode_usability(metadata.encode_usability(use)) == use
