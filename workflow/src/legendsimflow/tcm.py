# Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations

import awkward as ak
import numpy as np


def merge_stp_n_opt_tcms(tcm_stp, tcm_opt, *, scintillator_uid):
    """Merge `tcm_opt` rows into `tcm_stp` at the scintillator uid.

    For each `axis=0` row of `tcm_stp`, if `tcm_stp.table_key` contains
    `scintillator_uid`, replace that single element by splicing in the next row
    of `tcm_opt.table_key`. The same splice is applied to `row_in_table` using
    the corresponding `tcm_opt.row_in_table`, preserving alignment between
    `table_key[i][j]` and `row_in_table[i][j]`.

    Parameters
    ----------
    tcm_stp, tcm_opt
        Awkward record arrays with fields `table_key` and `row_in_table`.
    scintillator_uid
        Scalar value in `tcm_stp.table_key` marking where to splice in
        `tcm_opt`, i.e. the UID of the scintillator table.

    Returns
    -------
    ak.Array
        Record array with the same length as `tcm_stp`.
    """
    stp_k = tcm_stp.table_key
    stp_r = tcm_stp.row_in_table

    is_ph = stp_k == scintillator_uid
    n_ph_per_row = ak.sum(is_ph, axis=1)
    has_ph = n_ph_per_row > 0

    if ak.any(n_ph_per_row > 1):
        msg = "found multiple scintillator_uid placeholders in a single tcm_stp row"
        raise ValueError(msg)

    n_rows_with_ph = int(ak.sum(has_ph))
    if len(tcm_opt) != n_rows_with_ph:
        msg = (
            f"len(tcm_opt) ({len(tcm_opt)}) must equal number of tcm_stp rows containing "
            f"scintillator_uid ({n_rows_with_ph})"
        )
        raise ValueError(msg)

    # index of the placeholder in each row; None if absent
    pos_ph = ak.firsts(ak.local_index(stp_k)[is_ph])

    # map k-th row containing a placeholder -> tcm_opt[k]
    has_ph_np = ak.to_numpy(has_ph)
    opt_idx_np = np.cumsum(has_ph_np, dtype=np.int64) - 1
    opt_idx = ak.mask(ak.Array(opt_idx_np), has_ph)

    empty = ak.Array([[]])
    opt_k = ak.where(has_ph, tcm_opt.table_key[opt_idx], empty)
    opt_r = ak.where(has_ph, tcm_opt.row_in_table[opt_idx], empty)

    # splice per-row
    idx = ak.local_index(stp_k)
    before = idx < pos_ph
    after = idx > pos_ph

    left_k = ak.where(has_ph, stp_k[before], stp_k)
    right_k = ak.where(has_ph, stp_k[after], empty)
    merged_k = ak.where(has_ph, ak.concatenate([left_k, opt_k, right_k], axis=1), stp_k)

    left_r = ak.where(has_ph, stp_r[before], stp_r)
    right_r = ak.where(has_ph, stp_r[after], empty)
    merged_r = ak.where(has_ph, ak.concatenate([left_r, opt_r, right_r], axis=1), stp_r)

    return ak.zip({"table_key": merged_k, "row_in_table": merged_r}, depth_limit=1)
