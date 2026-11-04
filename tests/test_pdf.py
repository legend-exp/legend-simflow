from __future__ import annotations

import itertools

import hist
import lh5
import numpy as np
import pytest
from lgdo import Histogram

from legendsimflow.pdf import (
    build_pdf_grid,
    format_pdf_struct,
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
