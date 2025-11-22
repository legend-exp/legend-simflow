from __future__ import annotations

from pathlib import Path

from dbetto import AttrsDict

from legendsimflow import aggregate as agg


def test_simid_aggregates(fresh_config):
    config = fresh_config
    assert agg.get_simid_njobs(config, "stp", "birds_nest_K40") == 2
    assert isinstance(
        agg.gen_list_of_simid_inputs(config, "stp", "birds_nest_K40"),
        list,
    )
    assert isinstance(
        agg.gen_list_of_simid_outputs(config, "stp", "birds_nest_K40"),
        list,
    )

    config.benchmark.enabled = True
    config.benchmark.n_primaries.stp = 999

    assert agg.get_simid_njobs(config, "stp", "birds_nest_K40") == 1


def test_simid_harvesting(config):
    simids = agg.gen_list_of_all_simids(config, "stp")
    assert isinstance(simids, type({}.keys()))
    assert all(isinstance(s, str) for s in simids)


def test_simid_outputs(config):
    outputs = agg.gen_list_of_all_simid_outputs(config, "stp")
    assert isinstance(outputs, list)
    assert all(isinstance(s, Path) for s in outputs)
    assert len(outputs) == sum(
        [
            agg.get_simid_njobs(config, "stp", s)
            for s in agg.gen_list_of_all_simids(config, "stp")
        ]
    )


def test_process_simlist(config):
    targets = agg.process_simlist(
        config,
        simlist=["stp.birds_nest_K40", "stp.pen_plates_Ra224_to_Pb208"],
    )
    assert targets == agg.gen_list_of_simid_outputs(
        config, "stp", "birds_nest_K40"
    ) + agg.gen_list_of_simid_outputs(config, "stp", "pen_plates_Ra224_to_Pb208")


def test_dtmap_stuff(config):
    cry = agg.crystal_meta(
        config, config.metadata.hardware.detectors.germanium.diodes.V99000A
    )
    assert isinstance(cry, AttrsDict)
    assert cry.name == "000"
    assert cry.order == "99"

    assert agg.start_key(config, "l200-p02-r005-phy") == "20220602T000000Z"

    assert agg.gen_list_of_hpges_valid_for_dtmap(config, "l200-p02-r005-phy") == [
        "V99000A"
    ]
