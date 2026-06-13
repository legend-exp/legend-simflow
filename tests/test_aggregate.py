from __future__ import annotations

import logging
from pathlib import Path

import pytest
from dbetto import AttrsDict

import legendsimflow.aggregate as agg_mod
from legendsimflow import aggregate as agg
from legendsimflow.exceptions import SimflowConfigError
from legendsimflow.metadata import get_tier_settings


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
    assert len(simids) == 11


def test_simid_outputs(config):
    outputs = agg.gen_list_of_all_simid_outputs(config, "stp")
    assert isinstance(outputs, list)
    assert all(isinstance(s, Path) for s in outputs)
    assert len(outputs) == sum(
        [agg.get_simid_njobs(config, s) for s in agg.gen_list_of_all_simids(config)]
    )


def test_process_simlist_is_cumulative(config):
    # evt must include vtx/stp/opt/hit outputs for the same simid
    simid = "birds_nest_K40"
    targets_evt = agg.process_simlist(config, simlist=[f"evt.{simid}"])

    targets_stp = agg.process_simlist(config, simlist=[f"stp.{simid}"])
    targets_opt = agg.process_simlist(config, simlist=[f"opt.{simid}"])
    targets_hit = agg.process_simlist(config, simlist=[f"hit.{simid}"])

    assert targets_evt != []
    assert set(targets_stp).issubset(targets_evt)
    assert set(targets_opt).issubset(targets_evt)
    assert set(targets_hit).issubset(targets_evt)

    assert not set(targets_evt).issubset(targets_hit)
    assert not set(targets_evt).issubset(targets_opt)

    # cvt must include evt outputs (and therefore also lower tiers)
    simid2 = "pen_plates_Ra224_to_Pb208"
    targets_cvt = agg.process_simlist(config, simlist=[f"cvt.{simid2}"])
    targets_evt2 = agg.process_simlist(config, simlist=[f"evt.{simid2}"])

    assert targets_cvt != []
    assert set(targets_evt2).issubset(targets_cvt)


def test_process_simlist_is_cumulative_make_steps(config):
    make_steps = ["stp", "evt"]
    simid = "birds_nest_K40"
    targets_evt = agg.process_simlist(
        config, simlist=[f"evt.{simid}"], make_steps=make_steps
    )
    targets_stp = agg.process_simlist(
        config, simlist=[f"stp.{simid}"], make_steps=["stp"]
    )
    targets_opt = agg.process_simlist(
        config, simlist=[f"opt.{simid}"], make_steps=["opt"]
    )
    targets_hit = agg.process_simlist(
        config, simlist=[f"hit.{simid}"], make_steps=["hit"]
    )

    assert targets_evt != []
    assert set(targets_stp).issubset(targets_evt)

    assert not set(targets_opt).issubset(targets_evt)
    assert not set(targets_hit).issubset(targets_evt)


@pytest.mark.parametrize("step", ["vtx", "par"])
def test_process_simlist_rejects_non_simid_steps(config, step):
    with pytest.raises(SimflowConfigError):
        agg.process_simlist(config, simlist=[f"{step}.birds_nest_K40"])


def test_hpge_harvesting(config):
    cry = agg.crystal_meta(
        config, config.metadata.hardware.detectors.germanium.diodes.V99000A
    )
    assert isinstance(cry, AttrsDict)
    assert cry.name == "000"
    assert cry.order == "99"

    assert agg.start_key(config, "l200-p02-r005-phy") == "20220602T000000Z"

    # r000-r002 have an empty skip list, so V99000A is present; r003+ have V99000A skip-listed
    assert agg.gen_list_of_hpges_valid_for_modeling(config, "l200-p02-r001-phy") == [
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
    # now returns dict[str, int] (hpge -> voltage)
    assert hpges["l200-p02-r000-phy"] == {"V99000A": 4200}


def test_runlist_harvesting(config):
    assert agg.gen_list_of_all_runids(config) == {
        f"l200-p02-r00{i}-phy" for i in range(8)
    }


def test_dtmap_stuff(config):
    runid = "l200-p02-r000-phy"
    simid = "stp.pen_plates_Ra224_to_Pb208"

    dtmaps = agg.gen_list_of_dtmaps(config, runid)
    assert len(dtmaps) == 1
    # check that the dtmap filename contains the voltage
    assert "4200V" in str(dtmaps[0])

    assert len(agg.gen_list_of_merged_dtmaps(config, simid)) == 1
    assert len(agg.gen_list_of_dtmap_plots_outputs(config, simid)) == 1

    # dtmap plots are produced as par-tier plots, so they flow through the
    # generic plots aggregator
    assert agg.gen_list_of_plots_outputs(
        config, "par", simid
    ) == agg.gen_list_of_dtmap_plots_outputs(config, simid)
    assert len(agg.gen_list_of_all_plots_outputs(config, "par")) >= 1


def test_par_plots_psd_gate(fresh_config):
    config = fresh_config
    simid = "pen_plates_Ra224_to_Pb208"

    # with PSD enabled (the default) the par tier yields the dtmap plots
    assert len(agg.gen_list_of_plots_outputs(config, "par", simid)) >= 1

    # disabling PSD in the hit tier drops them
    get_tier_settings(config, "hit")["simulate_psd"] = False
    assert agg.gen_list_of_plots_outputs(config, "par", simid) == []


def test_hpge_voltage_functions(config):
    runid = "l200-p02-r000-phy"

    # test get_hpge_voltage
    voltage = agg.get_hpge_voltage(config, "V99000A", runid)
    assert voltage == 4200
    assert isinstance(voltage, int)


def test_currmod_stuff(config):
    runid = "l200-p02-r000-phy"
    simid = "stp.pen_plates_Ra224_to_Pb208"

    assert len(agg.gen_list_of_currmods(config, runid)) == 1
    assert len(agg.gen_list_of_merged_currmods(config, simid)) == 1
    assert len(agg.gen_list_of_elecmods(config, runid)) == 1
    assert len(agg.gen_list_of_merged_elecmods(config, simid)) == 1
    assert len(agg.gen_list_of_currmod_plots_outputs(config, simid)) == 1
    assert len(agg.gen_list_of_superpulses(config)) == 1


def test_psl_stuff(config):
    runid = "l200-p02-r000-phy"
    simid = "stp.pen_plates_Ra224_to_Pb208"

    ideal_psls = agg.gen_list_of_ideal_psls(config, runid)
    assert len(ideal_psls) == 1
    assert "4200V" in str(ideal_psls[0])

    realistic_psls = agg.gen_list_of_realistic_psls(config, runid)
    assert len(realistic_psls) == 1
    assert runid in str(realistic_psls[0])

    realistic_merged_psls = agg.gen_list_of_merged_realistic_psls(config, simid)
    assert len(realistic_merged_psls) == 1
    assert runid in str(realistic_merged_psls[0])


def test_tier_evt_stuff(config):
    files = agg.gen_list_of_all_tier_cvt_outputs(config)
    assert len(files) == 11


def test_usability_harvesting(config):
    usability = agg.gen_list_of_all_usabilities(config)
    assert isinstance(usability, dict)
    for v in usability.values():
        assert isinstance(v, dict)

    assert sorted(usability.keys()) == [
        "l200-p02-r000-phy",
        "l200-p02-r001-phy",
        "l200-p02-r002-phy",
        "l200-p02-r003-phy",
        "l200-p02-r004-phy",
        "l200-p02-r005-phy",
        "l200-p02-r006-phy",
        "l200-p02-r007-phy",
    ]
    assert usability["l200-p02-r000-phy"] == {
        "V99000A": {"usability": "on", "psd_usability": "valid"},
        "B99000A": {"usability": "off", "psd_usability": "missing"},
    }


def test_psd_usability_unknown_value(config, monkeypatch, caplog):
    """Unknown psd.status.low_aoe values should warn and default to 'valid'."""
    original = agg_mod.encode_psd_usability

    # Make encode_psd_usability raise KeyError for any value except "valid"
    # (simulates an unexpected psd.status.low_aoe in the metadata)
    def patched_encode(val):
        if val != "valid":
            raise KeyError(val)
        return original(val)

    monkeypatch.setattr(agg_mod, "encode_psd_usability", patched_encode)

    with caplog.at_level(logging.WARNING, logger="legendsimflow.aggregate"):
        result = agg_mod.gen_list_of_all_usabilities(config)

    assert "unexpected psd.status.low_aoe" in caplog.text
    # B99000A has psd.status.low_aoe="missing" -> triggers warning -> defaults to "valid"
    assert result["l200-p02-r000-phy"]["B99000A"]["psd_usability"] == "valid"


# skip-list tests
# Fixtures: tests/dummyprod/inputs/simprod/config/pars/legend/geds/skip/
#   validity.yaml:
#     - valid_from 20220102T000000Z -> l200-p02-r000-T%-all-skip.yaml  (empty dict {})
#     - valid_from 20220402T000000Z -> l200-p02-r003-T%-all-skip.yaml  (V99000A, B99000A)


def test_skip_list_removes_detector(config, caplog):
    """A detector listed in the skip file for a run is absent from the output."""
    # r003 start_key 20220402T000000Z matches the r003-T%-all-skip.yaml that lists V99000A
    with caplog.at_level(logging.WARNING, logger="legendsimflow.aggregate"):
        hpges = agg.gen_list_of_hpges_valid_for_modeling(config, "l200-p02-r003-phy")
    assert "V99000A" not in hpges


def test_skip_list_time_validity_boundary(config):
    """A detector in the skip file for r003+ is included for r000-r002 but not r003+."""
    # r000-r002: validity points to the empty skip file
    for runid in ("l200-p02-r000-phy", "l200-p02-r001-phy", "l200-p02-r002-phy"):
        hpges = agg.gen_list_of_hpges_valid_for_modeling(config, runid)
        assert "V99000A" in hpges

    # r003-r007: validity points to the skip file listing V99000A
    for runid in (
        "l200-p02-r003-phy",
        "l200-p02-r004-phy",
        "l200-p02-r005-phy",
        "l200-p02-r006-phy",
        "l200-p02-r007-phy",
    ):
        hpges = agg.gen_list_of_hpges_valid_for_modeling(config, runid)
        assert "V99000A" not in hpges


def test_skip_list_empty_mapping_is_noop(config):
    """An empty skip mapping {} leaves the HPGe list unchanged.

    r000-r002 map to l200-p02-r000-T%-all-skip.yaml which is {}: V99000A
    must still appear in the output.
    """
    hpges = agg.gen_list_of_hpges_valid_for_modeling(config, "l200-p02-r000-phy")
    assert "V99000A" in hpges
    assert isinstance(hpges, list)


def test_skip_list_usability_gate_takes_precedence(config, caplog):
    """A detector excluded by the usability gate is absent even if also skip-listed.

    B99000A has usability='off' in the channelmap (excluded before the skip
    check) and is also listed in the r003-T%-all-skip file. It must remain absent
    and must NOT produce a skip-list WARNING (the usability gate fires first).
    """
    with caplog.at_level(logging.WARNING, logger="legendsimflow.aggregate"):
        hpges = agg.gen_list_of_hpges_valid_for_modeling(config, "l200-p02-r003-phy")

    assert "B99000A" not in hpges
    # The skip-list warning should NOT fire for B99000A since the usability gate
    # prunes it before the skip check runs
    assert not any(
        "B99000A" in r.message for r in caplog.records if r.levelname == "WARNING"
    )


def test_skip_list_all_hpges_reflects_exclusion(config):
    """gen_list_of_all_hpges_valid_for_modeling omits detectors for all skip-covered runs."""
    hpges = agg.gen_list_of_all_hpges_valid_for_modeling(config)
    # r000-r002 have V99000A (empty skip list); r003-r007 do not
    for runid in ("l200-p02-r000-phy", "l200-p02-r001-phy", "l200-p02-r002-phy"):
        assert "V99000A" in hpges[runid]
    for runid in (
        "l200-p02-r003-phy",
        "l200-p02-r004-phy",
        "l200-p02-r005-phy",
        "l200-p02-r006-phy",
        "l200-p02-r007-phy",
    ):
        assert "V99000A" not in hpges[runid]
