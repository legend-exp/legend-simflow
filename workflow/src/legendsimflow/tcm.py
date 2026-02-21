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

from collections.abc import Iterable
from pathlib import Path

import awkward as ak
import numpy as np
import pygama.evt
from lgdo import Table, lh5


def build_tcm(
    hit_files: str | Path | Iterable[str | Path], out_file: str | Path
) -> None:
    """Re-create the TCM table from remage.

    Use remage fields `evtid` and `t0` (the latter is assumed to be in
    nanoseconds) to build coincidences. The settings are identical to the
    remage built-in TCM settings.
    """
    if isinstance(hit_files, str | Path):
        hit_files = [hit_files]

    # use tables keyed by UID in the __by_uid__ group.  in this way, the
    # TCM will index tables by UID.  the coincidence criterium is based
    # on Geant4 event identifier and time of the hits
    # NOTE: uses the same time window as in build_hit() reshaping
    pygama.evt.build_tcm(
        [(str(f), r"hit/__by_uid__/*") for f in hit_files],  # input_tables
        ["evtid", "t0"],  # coin_cols
        hash_func=r"(?<=hit/__by_uid__/det)\d+",
        coin_windows=[0, 10_000],  # in ns
        out_file=str(out_file),
        wo_mode="write_safe",
    )


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

    return merge_stp_n_opt_tcms_chunk(
        tcm_stp, tcm_opt, scintillator_uid=scintillator_uid
    )


def merge_stp_n_opt_tcms_chunk(tcm_stp, tcm_opt, *, scintillator_uid):
    """Chunk-level implementation of :func:`merge_stp_n_opt_tcms`.

    This function assumes `tcm_opt` contains *exactly* as many rows as there are
    rows in `tcm_stp` that contain `scintillator_uid`, in the same order.
    """
    stp_k = tcm_stp.table_key
    stp_r = tcm_stp.row_in_table

    is_ph = stp_k == scintillator_uid
    n_ph_per_row = ak.sum(is_ph, axis=1)
    has_ph = n_ph_per_row > 0

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


def merge_stp_n_opt_tcms_to_lh5(
    stp_file: str | Path,
    opt_file: str | Path,
    out_file: str | Path,
    *,
    scintillator_uid: int,
    buffer_len: str | int = "50*MB",
) -> None:
    """Stream-merge STP and OPT TCMs and write unified TCM to disk in chunks.

    Iterates over `stp_file:/tcm` using :class:`lgdo.lh5.LH5Iterator`. For each
    chunk, reads only the required number of OPT TCM rows (those corresponding
    to STP rows containing the `scintillator_uid` placeholder) via `lh5.read_as`
    with explicit indices. The merged output is appended to `out_file:/tcm`.
    """
    opt_pos = 0
    out_wo_mode = "write_safe"

    it_stp = lh5.LH5Iterator(
        str(stp_file),
        "tcm",
        buffer_len=buffer_len,
    )

    for stp_chunk in it_stp:
        tcm_stp_chunk = stp_chunk.view_as("ak")

        stp_k = tcm_stp_chunk.table_key
        is_ph = stp_k == scintillator_uid
        n_ph_per_row = ak.sum(is_ph, axis=1)

        if ak.any(n_ph_per_row > 1):
            msg = "found multiple scintillator_uid placeholders in a single tcm_stp row"
            raise ValueError(msg)

        has_ph = n_ph_per_row > 0
        n_need = int(ak.sum(has_ph))

        if n_need > 0:
            opt_rows = np.arange(opt_pos, opt_pos + n_need, dtype=np.int64)
            tcm_opt_chunk = lh5.read_as(
                "tcm", str(opt_file), idx=opt_rows, library="ak"
            )
            if len(tcm_opt_chunk) != n_need:
                msg = (
                    "Number of rows read from tcm_opt does not match expected count: "
                    f"requested {n_need}, got {len(tcm_opt_chunk)}"
                )
                raise ValueError(msg)
            opt_pos += n_need
        else:
            tcm_opt_chunk = ak.zip(
                {"table_key": ak.Array([]), "row_in_table": ak.Array([])},
                depth_limit=1,
            )

        merged = merge_stp_n_opt_tcms_chunk(
            tcm_stp_chunk,
            tcm_opt_chunk,
            scintillator_uid=scintillator_uid,
        )

        lh5.write(Table(merged), "tcm", str(out_file), wo_mode=out_wo_mode)
        out_wo_mode = "append"

    n_opt = lh5.read_n_rows("tcm", str(opt_file))
    if opt_pos != n_opt:
        msg = f"len(tcm_opt) ({n_opt}) must equal number of tcm_stp rows containing scintillator_uid ({opt_pos})"
        raise ValueError(msg)
