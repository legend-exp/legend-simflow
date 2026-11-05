from __future__ import annotations

import awkward as ak
import numpy as np
import pytest

from legendsimflow.exceptions import SimflowConfigError
from legendsimflow.selections import (
    aoe_categories,
    aoe_flags,
    aoe_model,
    attach_modelable_flag,
    category_levels,
    fccd_by_uid,
    fccd_categories,
    has_hpge_multiplicity,
    is_clean_event,
    is_hpge_hit_above_thr,
    is_spms_coincident,
    lar_categories,
    one_hit,
)


def _assert_partition(categories):
    np.testing.assert_array_equal(np.sum(list(categories.values()), axis=0), 1)


def test_hits_and_multiplicity():
    # a 10 keV hit is below the 25 keV threshold, in an AC detector
    evt = ak.Array(
        {
            "geds": {
                "energy": [
                    [1000.0],
                    [800.0, 10.0],
                    [400.0, 300.0],
                    [20.0],
                    [1.0, 2.0, 3.0],
                ],
                "quality": {
                    "is_good_channel": [
                        [True],
                        [True, False],
                        [True, False],
                        [True],
                        [True, True, True],
                    ]
                },
            }
        }
    )
    is_hit = is_hpge_hit_above_thr(evt)
    assert is_hit.tolist() == is_hpge_hit_above_thr(evt, 25).tolist()
    assert is_hit.tolist() == [
        [True],
        [True, False],
        [True, True],
        [False],
        [False] * 3,
    ]

    # event 2 has a hit in an AC detector: never selected
    assert has_hpge_multiplicity(evt, 1).tolist() == [True, True, False, False, False]
    assert has_hpge_multiplicity(evt, 2).tolist() == [False] * 5
    assert has_hpge_multiplicity(evt, 0, 1).tolist() == [True, True, False, True, True]
    assert has_hpge_multiplicity(evt, 1, np.inf).tolist() == [
        True,
        True,
        False,
        False,
        False,
    ]
    # a higher threshold drops the 300 keV hit of event 2, the one in the AC detector
    assert has_hpge_multiplicity(evt, 1, energy_thr_kev=350).tolist() == [
        True,
        True,
        True,
        False,
        False,
    ]

    m1 = evt[has_hpge_multiplicity(evt, 1)]
    assert one_hit(m1).energy.tolist() == [1000.0, 800.0]


def test_is_clean_event():
    evt = ak.Array(
        {
            "trigger": {"is_forced": [False, True, False, False, False, False]},
            "coincident": {
                "puls": [False, False, True, False, False, False],
                "muon": [False, False, False, True, False, False],
                "muon_offline": [False, False, False, False, True, False],
            },
            "geds": {"quality": {"is_bb_like": [True, True, True, True, True, False]}},
        }
    )
    assert is_clean_event(evt).tolist() == [True] + [False] * 5


def test_aoe_model():
    assert (
        aoe_model({"simulate_psd": True, "simulate_psd_with_psl": True}) == "pulse_lib"
    )
    assert (
        aoe_model({"simulate_psd": True, "simulate_psd_with_psl": False})
        == "single_temp"
    )
    assert aoe_model({"simulate_psd": False, "simulate_psd_with_psl": False}) is None


def test_category_levels():
    levels = category_levels(True, "pulse_lib", 10)
    assert [list(level) for level in levels] == [
        ["is_spms_coincident", "not_spms_coincident"],
        ["is_high_aoe", "is_ss_aoe", "is_low_aoe", "is_invalid_aoe"],
        ["is_fccd_tagged", "not_fccd_tagged"],
    ]
    assert [list(level) for level in category_levels(False, "single_temp", None)] == [
        ["is_ss_aoe", "is_low_aoe", "is_invalid_aoe"]
    ]
    assert category_levels(False, None, None) == []


def test_lar_veto():
    # - event 0: below both thresholds
    # - event 1: the photoelectron sum reaches the threshold
    # - event 2: the multiplicity reaches the threshold
    evt = ak.Array({"spms": {"multiplicity": [2, 2, 4], "energy_sum": [1.5, 4.0, 1.0]}})
    flag = is_spms_coincident(evt)
    assert flag.tolist() == [False, True, True]
    assert is_spms_coincident(evt, 5, 5).tolist() == [False, False, False]

    categories = lar_categories(flag)
    assert categories["is_spms_coincident"].tolist() == [False, True, True]
    assert categories["not_spms_coincident"].tolist() == [True, False, False]


def _per_hit(values):
    # one hit per event; the second hit of event 0 is below threshold
    return [[values[0], False], *([v] for v in values[1:])]


# - event 0: single-site
# - event 1: multi-site
# - event 2: above the high-side cut
# - event 3: no simulated A/E (not modeled in data)
# - event 4: no valid PSD
_ENERGY = [[1000.0, 10.0], [1000.0], [1000.0], [1000.0], [1000.0]]
_IS_GOOD = _per_hit([True, True, True, True, False])
_IS_SINGLE_SITE = _per_hit([True, False, True, False, True])
_IS_HIGH = _per_hit([False, False, True, False, False])


def _sim_events():
    model = {
        "has_aoe": _per_hit([True, True, True, False, True]),
        "is_single_site": _IS_SINGLE_SITE,
        "is_high_aoe": _IS_HIGH,
    }
    return ak.Array(
        {
            "geds": {
                "energy": _ENERGY,
                "psd": {"is_good": _IS_GOOD, "pulse_lib": model, "single_temp": model},
            }
        }
    )


def _data_events():
    return ak.Array(
        {
            "trigger": {"period": [2] * 5, "run": [0] * 5},
            "geds": {
                "energy": _ENERGY,
                "rawid": [[1, 0], [1], [1], [2], [1]],
                "psd": {
                    "is_good": _IS_GOOD,
                    "low_aoe": {"is_single_site": _IS_SINGLE_SITE},
                    "high_aoe": {"is_bulk": [[not v for v in h] for h in _IS_HIGH]},
                },
            },
        }
    )


@pytest.mark.parametrize(
    "events",
    [
        _sim_events,
        lambda: ak.with_field(
            _data_events(),
            [[True, False], [True], [True], [False], [True]],
            ["geds", "psd", "is_modeled"],
        ),
    ],
    ids=["sim", "data"],
)
def test_aoe_categories(events):
    categories = aoe_categories(*aoe_flags(events(), "pulse_lib"))

    assert categories["is_ss_aoe"].tolist() == [True, False, False, False, False]
    assert categories["is_low_aoe"].tolist() == [False, True, False, False, False]
    assert categories["is_high_aoe"].tolist() == [False, False, True, False, False]
    assert categories["is_invalid_aoe"].tolist() == [False, False, False, True, True]
    _assert_partition(categories)

    # no high-side cut: event 2 is single-site
    categories = aoe_categories(*aoe_flags(events(), "single_temp"))
    assert "is_high_aoe" not in categories
    assert categories["is_ss_aoe"].tolist() == [True, False, True, False, False]
    _assert_partition(categories)


def test_aoe_flags_need_modelable_flag():
    with pytest.raises(ValueError, match="attach_modelable_flag"):
        aoe_flags(_data_events(), "pulse_lib")


def test_attach_modelable_flag(config):
    # the rawids of V99000A and B99000A in the p02 channel map; event 3 is in
    # run 1, where B99000A is modeled instead
    evt = _data_events()
    evt["trigger", "run"] = [0, 0, 0, 1, 0]
    evt["geds", "rawid"] = [
        [1234568, 1234569],
        [1234568],
        [1234568],
        [1234569],
        [1234569],
    ]
    is_modelable = {
        "l200-p02-r000-phy": {"V99000A": True, "B99000A": False},
        "l200-p02-r001-phy": {"V99000A": False, "B99000A": True},
    }
    evt = attach_modelable_flag(evt, config, is_modelable)

    assert evt.geds.psd.is_modeled.tolist() == [
        [True, False],
        [True],
        [True],
        [True],
        [False],
    ]

    with pytest.raises(SimflowConfigError, match=r"period 2 run 1 .* need one run"):
        attach_modelable_flag(evt, config, {"l200-p02-r000-phy": {}})


def test_fccd_categories():
    # FCCD 1 mm; the energy in the FCCD layer is 5, 15 and 0 keV
    evt = ak.Array(
        {
            "geds": {
                "energy": [[1000.0], [1000.0], [1000.0]],
                "surface": {
                    "edep": [[[5.0, 20.0]], [[10.0, 5.0]], [[50.0]]],
                    "distance_to_nplus": [[[0.5, 2.0]], [[0.1, 0.9]], [[1.5]]],
                },
            }
        }
    )
    categories = fccd_categories(evt, np.array([1.0, 1.0, 1.0]), fccd_thr_kev=10)

    assert categories["is_fccd_tagged"].tolist() == [False, True, False]
    _assert_partition(categories)


def test_fccd_by_uid(config, caplog):
    uids, fccd_mm = fccd_by_uid(
        config.metadata, {"V99000A": 5, "B99000A": 2}, surface_depth=3
    )

    assert uids.tolist() == [2, 5]
    assert fccd_mm.tolist() == pytest.approx([0.75, 0.75])

    with caplog.at_level("WARNING"):
        _, fccd_mm = fccd_by_uid(config.metadata, {"V99000A": 5}, surface_depth=0.5)
    assert fccd_mm.tolist() == [0.5]
    assert "exceeds surface_depth_in_mm" in caplog.text
