from __future__ import annotations

from lgdo import lh5

from legendsimflow import reboost as rutils


def test_remage_hit_range(legend_testdata):
    f_stp = legend_testdata["remage/th228-full-optional-v0_13.lh5"]

    n_rows = lh5.read_n_rows("stp/det1", f_stp)
    assert rutils.get_remage_hit_range(f_stp, "det1", 11, [0, -1]) == (0, n_rows)

    tcm = lh5.read_as("tcm", f_stp, "ak")
    assert rutils.get_remage_hit_range(f_stp, "det1", 11, [0, len(tcm)]) == (0, n_rows)
    assert rutils.get_remage_hit_range(f_stp, "det1", 11, [22, 98]) == (0, 3)
    assert rutils.get_remage_hit_range(f_stp, "det1", 11, [22, 97]) == (0, 2)
