from __future__ import annotations

import awkward as ak

from legendsimflow.tcm import merge_stp_n_opt_tcms


def test_merge_tcm_basic():
    tcm_stp = ak.zip(
        {
            "table_key": ak.Array([[-1, 1086400], [1113603, -1], [42]]),
            "row_in_table": ak.Array([[1, 0], [0, 2], [7]]),
        },
        depth_limit=1,
    )

    tcm_opt = ak.zip(
        {
            "table_key": ak.Array([[1052802, 1052803], [1054401]]),
            "row_in_table": ak.Array([[1, 1], [2]]),
        },
        depth_limit=1,
    )

    out = merge_stp_n_opt_tcms(tcm_stp, tcm_opt, scintillator_uid=-1)

    assert out.table_key[0].to_list() == [1052802, 1052803, 1086400]
    assert out.row_in_table[0].to_list() == [1, 1, 0]

    assert out.table_key[1].to_list() == [1113603, 1054401]
    assert out.row_in_table[1].to_list() == [0, 2]

    assert out.table_key[2].to_list() == [42]
    assert out.row_in_table[2].to_list() == [7]

    assert len(out) == len(tcm_stp)
    assert not ak.any(out.table_key == -1)
