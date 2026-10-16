from __future__ import annotations

import hist
import lh5
import numpy as np
import pytest

from legendsimflow.pdf_grid import build_pdf_grid


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


def test_build_pdf_grid_uneven_axis():
    with pytest.raises(ValueError, match="does not go"):
        build_pdf_grid(_pdf, {"a": (0, 1, 0.3)})


def test_build_pdf_grid_different_binning():
    def pdf(pars):
        return hist.new.Reg(10 + int(pars["a"]), 0, 40).Double()

    with pytest.raises(ValueError, match="different binning"):
        build_pdf_grid(pdf, {"a": (0, 1, 1)})


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


def test_build_pdf_grid_changing_keys():
    def pdf(pars):
        return {f"pdf_{int(pars['a'])}": _pdf(pars | {"b": 0})}

    with pytest.raises(ValueError, match="differ from those"):
        build_pdf_grid(pdf, {"a": (0, 1, 1)})
