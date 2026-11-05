from __future__ import annotations

import awkward as ak
import numpy as np
import pytest

from legendsimflow.selections import (
    aoe_categories,
    aoe_model,
    category_levels,
    data_aoe_flags,
    data_lar_veto,
    fccd_by_uid,
    fccd_categories,
    is_clean_data_event,
    lar_categories,
    select_hits,
    sim_aoe_flags,
)


def _assert_partition(categories):
    np.testing.assert_array_equal(np.sum(list(categories.values()), axis=0), 1)


def test_select_hits():
    # a 10 keV hit is below the 25 keV threshold, in an AC detector
    evt = ak.Array(
        {
            "geds": {
                "energy": [[1000.0], [800.0, 10.0], [300.0, 400.0], [20.0]],
                "quality": {
                    "is_good_channel": [[True], [True, False], [True, False], [True]]
                },
            }
        }
    )
    is_hit, multiplicity, is_good = select_hits(evt, 25)

    assert is_hit.tolist() == [[True], [True, False], [True, True], [False]]
    assert multiplicity.tolist() == [1, 1, 2, 0]
    assert is_good.tolist() == [True, True, False, True]


def test_is_clean_data_event():
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
    assert is_clean_data_event(evt).tolist() == [True] + [False] * 5


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


def test_lar_categories():
    categories = lar_categories(np.array([True, False, False]))

    assert categories["is_spms_coincident"].tolist() == [True, False, False]
    assert categories["not_spms_coincident"].tolist() == [False, True, True]


def test_data_lar_veto():
    # - event 0: below both thresholds
    # - event 1: the photoelectron sum reaches the threshold
    # - event 2: the multiplicity reaches the threshold
    evt = ak.Array({"spms": {"multiplicity": [2, 2, 4], "energy_sum": [1.5, 4.0, 1.0]}})
    flag = data_lar_veto(evt, multiplicity_thr=4, energy_sum_thr_pe=4)
    assert flag.tolist() == [False, True, True]


def _per_hit(values):
    # one hit per event; the second hit of event 0 is not selected
    return [[values[0], False], *([v] for v in values[1:])]


# - event 0: single-site
# - event 1: multi-site
# - event 2: above the high-side cut
# - event 3: no simulated A/E (not modeled in data)
# - event 4: no valid PSD
_IS_HIT = ak.Array([[True, False], [True], [True], [True], [True]])
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
                "psd": {"is_good": _IS_GOOD, "pulse_lib": model, "single_temp": model}
            }
        }
    )


def _data_events():
    return ak.Array(
        {
            "geds": {
                "rawid": [[1, 0], [1], [1], [2], [1]],
                "psd": {
                    "is_good": _IS_GOOD,
                    "low_aoe": {"is_single_site": _IS_SINGLE_SITE},
                    "high_aoe": {"is_bulk": [[not v for v in h] for h in _IS_HIGH]},
                },
            }
        }
    )


@pytest.mark.parametrize(
    "flags",
    [
        lambda model: sim_aoe_flags(_sim_events(), _IS_HIT, model),
        lambda model: data_aoe_flags(_data_events(), _IS_HIT, model, [1]),
    ],
    ids=["sim", "data"],
)
def test_aoe_categories(flags):
    categories = aoe_categories(*flags("pulse_lib"))

    assert categories["is_ss_aoe"].tolist() == [True, False, False, False, False]
    assert categories["is_low_aoe"].tolist() == [False, True, False, False, False]
    assert categories["is_high_aoe"].tolist() == [False, False, True, False, False]
    assert categories["is_invalid_aoe"].tolist() == [False, False, False, True, True]
    _assert_partition(categories)

    # no high-side cut: event 2 is single-site
    categories = aoe_categories(*flags("single_temp"))
    assert "is_high_aoe" not in categories
    assert categories["is_ss_aoe"].tolist() == [True, False, True, False, False]
    _assert_partition(categories)


def test_fccd_categories():
    # FCCD 1 mm; the energy in the FCCD layer is 5, 15 and 0 keV
    evt = ak.Array(
        {
            "geds": {
                "surface": {
                    "edep": [[[5.0, 20.0]], [[10.0, 5.0]], [[50.0]]],
                    "distance_to_nplus": [[[0.5, 2.0]], [[0.1, 0.9]], [[1.5]]],
                }
            }
        }
    )
    is_hit = ak.Array([[True], [True], [True]])
    categories = fccd_categories(
        evt, is_hit, np.array([1.0, 1.0, 1.0]), energy_thr_kev=10
    )

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
