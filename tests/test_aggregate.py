from __future__ import annotations

from pathlib import Path

from dbetto import AttrsDict

from legendsimflow import aggregate as agg


def test_simid_aggregates(fresh_config):
    config = fresh_config
    assert agg.get_simid_njobs(config, "birds_nest_K40") == 2

    val = agg.gen_list_of_simid_inputs(config, "stp", "birds_nest_K40")
    assert isinstance(val, list)
    assert len(val) == 1

    val = agg.gen_list_of_simid_outputs(config, "stp", "birds_nest_K40")
    assert isinstance(val, list)
    assert len(val) == 2

    val = agg.gen_list_of_simid_outputs(config, "vtx", "birds_nest_K40")
    assert isinstance(val, list)
    assert len(val) == 2

    config.benchmark.enabled = True
    config.benchmark.n_primaries.stp = 999

    assert agg.get_simid_njobs(config, "birds_nest_K40") == 1


def test_simid_harvesting(config):
    simids = agg.gen_list_of_all_simids(config)
    assert isinstance(simids, type({}.keys()))
    assert all(isinstance(s, str) for s in simids)
    assert len(simids) == 9


def test_simid_outputs(config):
    outputs = agg.gen_list_of_all_simid_outputs(config, "stp")
    assert isinstance(outputs, list)
    assert all(isinstance(s, Path) for s in outputs)
    assert len(outputs) == sum(
        [agg.get_simid_njobs(config, s) for s in agg.gen_list_of_all_simids(config)]
    )


def test_process_simlist(config):
    for tier in ("vtx", "stp", "opt", "hit", "evt"):
        targets = agg.process_simlist(
            config,
            simlist=[f"{tier}.birds_nest_K40", f"{tier}.pen_plates_Ra224_to_Pb208"],
        )
        assert targets != []


def test_hpge_harvesting(config):
    cry = agg.crystal_meta(
        config, config.metadata.hardware.detectors.germanium.diodes.V99000A
    )
    assert isinstance(cry, AttrsDict)
    assert cry.name == "000"
    assert cry.order == "99"

    assert agg.start_key(config, "l200-p02-r005-phy") == "20220602T000000Z"

    assert agg.gen_list_of_hpges_valid_for_modeling(config, "l200-p02-r005-phy") == [
        "V99000A"
    ]

    hpges = agg.gen_list_of_all_hpges_valid_for_modeling(config)
    assert isinstance(hpges, dict)
    assert sorted(hpges.keys()) == [
        "l200-p02-r000-phy",
        "l200-p02-r001-phy",
        "l200-p02-r002-phy",
        "l200-p02-r003-phy",
        "l200-p02-r004-phy",
        "l200-p02-r005-phy",
        "l200-p02-r006-phy",
        "l200-p02-r007-phy",
    ]
    assert hpges["l200-p02-r000-phy"] == ["V99000A"]


def test_runlist_harvesting(config):
    assert agg.gen_list_of_all_runids(config) == {
        f"l200-p02-r00{i}-phy" for i in range(8)
    }


def test_dtmap_stuff(config):
    runid = "l200-p02-r000-phy"
    simid = "stp.pen_plates_Ra224_to_Pb208"

    assert len(agg.gen_list_of_dtmaps(config, runid)) == 1
    assert len(agg.gen_list_of_merged_dtmaps(config, simid)) == 1
    assert len(agg.gen_list_of_dtmap_plots_outputs(config, simid)) == 1

    plots = agg.gen_list_of_all_dtmap_plots_outputs(config)
    assert isinstance(plots, set)
    assert len(plots) == 8


def test_currmod_stuff(config):
    runid = "l200-p02-r000-phy"
    simid = "stp.pen_plates_Ra224_to_Pb208"

    assert len(agg.gen_list_of_currmods(config, runid)) == 1
    assert len(agg.gen_list_of_merged_currmods(config, simid)) == 1
    assert len(agg.gen_list_of_currmod_plots_outputs(config, simid)) == 1

    plots = agg.gen_list_of_all_currmod_plots_outputs(config)
    assert isinstance(plots, set)
    assert len(plots) == 8


def test_tier_evt_stuff(config):
    files = agg.gen_list_of_all_concat_tier_evt_outputs(config)
    assert len(files) == 9
