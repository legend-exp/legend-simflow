from __future__ import annotations

import itertools

import awkward as ak
import hist
import lh5
import numpy as np
import pytest
from lgdo import Histogram
from reboost.math.functions import piecewise_linear_activeness as reference

from legendsimflow.pdf import (
    aoe_categories,
    build_pdf_grid,
    collected_energy,
    fccd_by_uid,
    fccd_categories,
    flatten_surface_steps,
    format_pdf_struct,
    lar_categories,
    piecewise_linear_activeness,
    read_pdf,
    validate_grid,
)


def _pdf(pars):
    # one count at x = a + 10 b, so each grid point is recognisable
    h = hist.new.Reg(40, 0, 40).Double()
    h.fill([pars["a"] + 10 * pars["b"]])
    return h


def test_build_pdf_grid(tmp_path):
    axes = {"a": (0, 2, 1), "b": (0, 1, 0.5)}
    grid = build_pdf_grid(_pdf, axes)

    info = grid["grid_info"]
    assert info["axis_0"]["name"].value == "a"
    assert info["axis_1"]["name"].value == "b"
    assert info["axis_1"]["values"]["step"].value == 0.5

    assert sorted(grid["grid"].keys()) == ["a_0", "a_1", "a_2"]
    assert sorted(grid["grid"]["a_0"].keys()) == ["b_0", "b_1", "b_2"]

    # a = 2, b = 0.5 -> x = 7
    weights = grid["grid"]["a_2"]["b_1"].weights.nda
    assert weights[7] == 1
    assert weights.sum() == 1

    lh5.write(grid, "pdf_grid", tmp_path / "grid.lh5")
    read_back = lh5.read("pdf_grid/grid/a_2/b_1", tmp_path / "grid.lh5")
    assert np.array_equal(read_back.weights.nda, weights)


def test_build_pdf_grid_several_pdfs(tmp_path):
    def pdf(pars):
        return {"hit/all": _pdf(pars), "mul1/all": _pdf({"a": 0, "b": pars["b"]})}

    grids = build_pdf_grid(pdf, {"a": (0, 2, 1), "b": (0, 1, 0.5)})
    assert sorted(grids) == ["hit/all", "mul1/all"]
    assert grids["hit/all"]["grid"]["a_2"]["b_1"].weights.nda[7] == 1
    assert grids["mul1/all"]["grid"]["a_2"]["b_1"].weights.nda[5] == 1

    for path, grid in grids.items():
        lh5.write(grid, f"pdf_grid/{path}", tmp_path / "grid.lh5")
    assert (
        lh5.read("pdf_grid/mul1/all/grid/a_2/b_1", tmp_path / "grid.lh5").weights.nda[5]
        == 1
    )


def test_validate_grid(tmp_path):
    grid = build_pdf_grid(_pdf, {"a": (0, 2, 1), "b": (0, 1, 0.5)})
    lh5.write(grid, "pdf_grid", tmp_path / "grid.lh5")
    validate_grid(lh5.read("pdf_grid", tmp_path / "grid.lh5"))

    grid["grid"]["a_1"].remove_field("b_2")
    with pytest.raises(ValueError, match="a_1/b_2 has no histogram"):
        validate_grid(grid)


def test_validate_grid_uneven_axis():
    with pytest.raises(ValueError, match="does not go"):
        build_pdf_grid(_pdf, {"a": (0, 1, 0.3), "b": (0, 1, 1)})


def test_validate_grid_different_binning():
    def pdf(pars):
        return hist.new.Reg(10 + int(pars["a"]), 0, 40).Double()

    with pytest.raises(ValueError, match="different binning"):
        build_pdf_grid(pdf, {"a": (0, 1, 1)})


def _assert_partition(categories):
    np.testing.assert_array_equal(np.sum(list(categories.values()), axis=0), 1)


def test_lar_categories():
    evt = ak.Array({"coincident": {"spms": [True, False, False]}})
    categories = lar_categories(evt)

    assert categories["is_spms_coincident"].tolist() == [True, False, False]
    assert categories["not_spms_coincident"].tolist() == [False, True, True]


def _psd_hits():
    # - hit 0: single-site
    # - hit 1: multi-site
    # - hit 2: above the high-side cut
    # - hit 3: no simulated A/E
    # - hit 4: no valid PSD in data
    model = {
        "has_aoe": [True, True, True, False, True],
        "is_single_site": [True, False, True, False, True],
        "is_high_aoe": [False, False, True, False, False],
    }
    return ak.Array(
        {
            "is_good": [True, True, True, True, False],
            "pulse_lib": model,
            "single_temp": model,
        }
    )


def test_aoe_categories_pulse_lib():
    categories = aoe_categories(_psd_hits(), "pulse_lib")

    assert categories["is_ss_aoe"].tolist() == [True, False, False, False, False]
    assert categories["is_low_aoe"].tolist() == [False, True, False, False, False]
    assert categories["is_high_aoe"].tolist() == [False, False, True, False, False]
    assert categories["is_invalid_aoe"].tolist() == [False, False, False, True, True]
    _assert_partition(categories)


def test_aoe_categories_single_temp():
    categories = aoe_categories(_psd_hits(), "single_temp")

    # no high-side cut: event 2 is single-site
    assert "is_high_aoe" not in categories
    assert categories["is_ss_aoe"].tolist() == [True, False, True, False, False]
    _assert_partition(categories)


def _geds_surface():
    # the second hit of event 1 (UID 2) has no depositions near the surface
    return ak.Array(
        {
            "energy": [[1000.0], [800.0, 300.0]],
            "rawid": [[1], [1, 2]],
            "surface": {
                "edep": [[[10.0, 20.0]], [[5.0], []]],
                "distance_to_nplus": [[[0.5, 2.0]], [[0.1], []]],
                "energy_bulk": [[970.0], [795.0, 300.0]],
            },
        }
    )


def test_flatten_surface_steps():
    steps = flatten_surface_steps(_geds_surface(), np.array([1, 2]))

    assert steps["energy_bulk"].tolist() == [970.0, 795.0, 300.0]
    assert steps["hit"].tolist() == [0, 0, 1]
    assert steps["det"].tolist() == [0, 0, 0]
    assert steps["edep"].tolist() == [10.0, 20.0, 5.0]


def test_fccd_categories():
    # FCCD 1 mm; the energy in the FCCD layer is 10, 5 and 0 keV
    steps = flatten_surface_steps(_geds_surface(), np.array([1, 2]))
    categories = fccd_categories(steps, np.array([1.0, 1.0]), energy_thr_kev=6)

    assert categories["is_fccd_tagged"].tolist() == [True, False, False]
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


@pytest.mark.parametrize("dlf", [0, 0.5, 1])
def test_piecewise_linear_activeness(dlf):
    distance = np.array([0, 0.2, 0.5, 0.9, 1.0, 1.5])
    assert piecewise_linear_activeness(distance, 1.0, dlf) == pytest.approx(
        reference(ak.Array(distance), 1.0, dlf).to_numpy()
    )


def test_collected_energy():
    steps = flatten_surface_steps(_geds_surface(), np.array([1, 2]))

    # activeness at 0.5 mm and 0.1 mm, FCCD 1 mm: 0 and 0 with DLF 0.5, 0.5 and
    # 0.1 with DLF 0
    collected = collected_energy(steps, np.array([1.0, 1.0]), 0.5)
    assert collected.tolist() == pytest.approx([990.0, 795.0, 300.0])
    collected = collected_energy(steps, np.array([1.0, 1.0]), 0)
    assert collected.tolist() == pytest.approx([995.0, 795.5, 300.0])

    # FCCD 0.4 mm: the deposition at 0.5 mm is fully collected
    collected = collected_energy(steps, np.array([0.4, 1.0]), 0.5)
    assert collected.tolist() == pytest.approx([1000.0, 795.0, 300.0])


def test_format_pdf_struct():
    histograms = {
        "mul1/a/b": {"all": hist.new.Reg(10, 0, 10).Double()},
        "mul1/a/c": {"all": hist.new.Reg(10, 0, 10).Double()},
    }
    struct = format_pdf_struct(histograms, {"mul1/a/b": "B", "mul1/a/c": "C"})

    assert struct["mul1"]["a"]["b"]["all"].attrs["description"] == "B"
    assert struct["mul1"]["a"]["c"]["all"].attrs["description"] == "C"


def test_read_pdf(tmp_path):
    # one event at energy n keV in each category n, in both groups
    lar = ("is_spms_coincident", "not_spms_coincident")
    aoe = ("is_ss_aoe", "is_low_aoe")
    histograms, descriptions = {}, {}
    for n, (lar_cat, aoe_cat) in enumerate(itertools.product(lar, aoe)):
        path = f"mul1/{lar_cat}/{aoe_cat}"
        histograms[path] = {}
        for g in ("all", "icpc"):
            histograms[path][g] = hist.new.Reg(10, 0, 10).Double()
            histograms[path][g].fill([n])
        descriptions[path] = path
    struct = format_pdf_struct(histograms, descriptions)
    struct.add_field(
        "mul2/all", Histogram(hist.new.Reg(10, 0, 10).Reg(10, 0, 10).Double())
    )
    pdf_file = tmp_path / "pdf.lh5"
    lh5.write(struct, "pdf", pdf_file)

    def _counts(selection, **kwargs):
        return read_pdf(selection, pdf_file, **kwargs).values()[:4].tolist()

    assert _counts("mul1") == [1, 1, 1, 1]
    assert _counts("mul1/not_spms_coincident") == [0, 0, 1, 1]
    assert _counts("mul1/*/is_ss_aoe") == [1, 0, 1, 0]
    assert _counts(["mul1/*/is_ss_aoe", "mul1/is_spms_coincident"]) == [1, 1, 1, 0]
    assert _counts("mul1/not_spms_coincident/is_low_aoe", group="icpc") == [0, 0, 0, 1]

    with pytest.raises(ValueError, match="no histogram"):
        read_pdf("mul1/*/is_high_aoe", pdf_file)
