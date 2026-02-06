from __future__ import annotations

import awkward as ak
from lgdo import lh5

from legendsimflow import awkward


def test_ak_isin_basic():
    a = ak.Array([[1, 2, 3], [2, 4]])
    b = ak.Array([2, 3])

    result = awkward.ak_isin(a, b)
    expected = ak.Array([[False, True, True], [True, False]])

    assert ak.to_list(result) == ak.to_list(expected)


def test_structure_preserved():
    a = ak.Array([[], [1], [1, 2, 3]])
    b = ak.Array([2])

    result = awkward.ak_isin(a, b)
    assert ak.num(result).tolist() == ak.num(a).tolist()


def test_tcm(legend_testdata):
    tcm = lh5.read_as(
        "tcm",
        legend_testdata["remage/th228-full-optional-v0_13.lh5"],
        library="ak",
        n_rows=3,
    )
    key1 = awkward.ak_isin(tcm.table_key, [102, 1])
    assert key1.tolist() == [
        [False, True],
        [False, True],
        [True, False, False, True],
    ]
    key2 = awkward.ak_isin(tcm.table_key, [101, 2])
    assert key2.tolist() == [[True, False], [True, False], [False, True, True, False]]

    tcm1 = tcm[key1]
    tcm2 = tcm[key2]

    tcm_new = ak.sort(ak.concatenate([tcm1, tcm2], axis=-1), axis=-1)
    tcm_sort = ak.sort(tcm, axis=-1)

    assert tcm_new.table_key.tolist() == tcm_sort.table_key.tolist()
    assert tcm_new.row_in_table.tolist() == tcm_sort.row_in_table.tolist()
