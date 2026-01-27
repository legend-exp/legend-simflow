from __future__ import annotations

from lgdo import lh5

from legendsimflow import reboost as rutils


def test_remage_hit_range(legend_testdata):
    f_stp = legend_testdata["remage/th228-full-optional-v0_13.lh5"]

    tcm = lh5.read_as("tcm", f_stp, library="ak")

    for det, uid in zip(
        ["det1", "det2", "scint1", "scint2", "optdet1", "optdet2"],
        [11, 12, 1, 2, 101, 102],
        strict=True,
    ):
        n_rows = lh5.read_n_rows(f"stp/{det}", f_stp)

        assert rutils.get_remage_hit_range(tcm, det, uid, [0, -1]) == (0, n_rows)

        tcm = lh5.read_as("tcm", f_stp, "ak")

        assert rutils.get_remage_hit_range(tcm, det, uid, [0, len(tcm)]) == (0, n_rows)

    det = "det1"
    uid = 11
    assert rutils.get_remage_hit_range(tcm, "det1", 11, [22, 98]) == (0, 3)
    assert rutils.get_remage_hit_range(tcm, "det1", 11, [22, 97]) == (0, 2)
