from __future__ import annotations

import hist
import lh5
import numpy as np
import pytest
from lgdo import Scalar

from legendsimflow.pdf import (
    build_pdf_grid,
    fill_histograms,
    format_pdf_struct,
    group_detectors,
    merge_pdf_files,
    new_histograms,
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


def test_group_detectors(caplog):
    detectors = {"V01": 1, "B02": 2, "V03": 3}
    groups = group_detectors({"icpc": "V.*", "bege": "B.*"}, detectors)
    assert {g: ids.tolist() for g, ids in groups.items()} == {
        "icpc": [1, 3],
        "bege": [2],
        "all": [1, 2, 3],
    }
    assert group_detectors(None, detectors)["all"].tolist() == [1, 2, 3]

    with caplog.at_level("WARNING"):
        groups = group_detectors({"all": "V.*", "coax": "C.*"}, detectors)
    assert groups["all"].tolist() == [1, 2, 3]
    assert groups["coax"].tolist() == []
    assert "overridden by the implicit all-detector group" in caplog.text
    assert "'coax' matches no detectors" in caplog.text


_LEVELS = [{"a": "A", "b": "B"}, {"c": "C", "d": "D"}]


def test_fill_histograms():
    groups = {"all": np.array([1, 2]), "icpc": np.array([1])}
    histograms, mul2_hist = new_histograms(_LEVELS, groups)
    assert sorted(histograms) == [("a", "c"), ("a", "d"), ("b", "c"), ("b", "d")]
    assert mul2_hist.ndim == 2

    # events: (a, c) in 1, (a, d) in 2, (b, c) in 1
    energy = np.array([100.0, 200.0, 300.0])
    detector = np.array([1, 2, 1])
    in_category = {
        "a": np.array([True, True, False]),
        "b": np.array([False, False, True]),
        "c": np.array([True, False, True]),
        "d": np.array([False, True, False]),
    }
    fill_histograms(histograms, energy, detector, in_category, groups)

    assert histograms["a", "c"]["all"][100j] == 1
    assert histograms["a", "c"]["icpc"][100j] == 1
    assert histograms["a", "d"]["all"][200j] == 1
    assert histograms["a", "d"]["icpc"].sum() == 0
    assert histograms["b", "c"]["all"].sum() == 1
    assert histograms["b", "d"]["all"].sum() == 0


def test_format_pdf_struct():
    histograms, mul2_hist = new_histograms(_LEVELS, ["all"])
    struct = format_pdf_struct(histograms, _LEVELS, mul2_hist)

    assert (
        struct["mul1"]["a"]["c"]["all"].attrs["description"] == "multiplicity 1, A, C"
    )
    assert (
        struct["mul1"]["b"]["d"]["all"].attrs["description"] == "multiplicity 1, B, D"
    )
    assert "mul2" in struct
    assert "mul2" not in format_pdf_struct(histograms, _LEVELS)


def test_read_pdf(tmp_path):
    # one event at energy n keV in each category n, in both groups
    levels = [
        {"is_spms_coincident": "", "not_spms_coincident": ""},
        {"is_ss_aoe": "", "is_low_aoe": ""},
    ]
    histograms, mul2_hist = new_histograms(levels, ["all", "icpc"])
    for n, per_group in enumerate(histograms.values()):
        for h in per_group.values():
            h.fill([n])
    pdf_file = tmp_path / "pdf.lh5"
    lh5.write(format_pdf_struct(histograms, levels, mul2_hist), "pdf", pdf_file)

    def _counts(selection, **kwargs):
        return read_pdf(selection, pdf_file, **kwargs).values()[:4].tolist()

    assert _counts("mul1") == [1, 1, 1, 1]
    assert _counts("mul1/not_spms_coincident") == [0, 0, 1, 1]
    assert _counts("mul1/*/is_ss_aoe") == [1, 0, 1, 0]
    assert _counts(["mul1/*/is_ss_aoe", "mul1/is_spms_coincident"]) == [1, 1, 1, 0]
    assert _counts("mul1/not_spms_coincident/is_low_aoe", group="icpc") == [0, 0, 0, 1]

    with pytest.raises(ValueError, match="no histogram"):
        read_pdf("mul1/*/is_high_aoe", pdf_file)


def test_merge_pdf_files(tmp_path):
    files = []
    for n in (1, 2):
        histograms, mul2_hist = new_histograms(_LEVELS, ["all"])
        histograms["a", "c"]["all"].fill([10.0] * n)
        mul2_hist.fill([100.0], [200.0])
        files.append(tmp_path / f"{n}.lh5")
        lh5.write(format_pdf_struct(histograms, _LEVELS, mul2_hist), "pdf", files[-1])
        lh5.write(Scalar(1000.0 * n), "livetime_in_s", files[-1], wo_mode="append")

    merged = tmp_path / "merged.lh5"
    merge_pdf_files(files, merged)

    assert read_pdf("mul1/a/c", merged)[10j] == 3
    assert read_pdf("mul1", merged).sum() == 3
    assert lh5.read_as("pdf/mul2/all", merged, "hist").sum() == 2
    assert (
        lh5.read("pdf/mul1/a/c/all", merged).attrs["description"]
        == "multiplicity 1, A, C"
    )
    assert lh5.read("livetime_in_s", merged).value == 3000
    assert "nr_sim_events" not in lh5.ls(merged)
