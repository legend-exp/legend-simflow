# Copyright (C) 2026 Luigi Pertoldi <gipert@pm.me>
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

import runpy
from pathlib import Path
from types import SimpleNamespace

import matplotlib as mpl

mpl.use("Agg")

import lh5
import matplotlib.pyplot as plt
import numpy as np
import pytest
from dbetto import AttrsDict
from lgdo import Array, Scalar, Struct, Table, VectorOfVectors
from snakemake.iocontainers import InputFiles

import legendsimflow

# the cvt-observables plot is a pure Snakemake script (no CLI ``main()``); it is
# driven here via ``runpy`` with a mock ``snakemake`` global, mirroring how
# Snakemake injects it at run time.
_PLOT_SCRIPT = (
    Path(legendsimflow.__file__).parent
    / "scripts"
    / "plots"
    / "tier_cvt_observables.py"
)


def _write_cvt_root(path: Path, names: list[str], uids: list[int]) -> None:
    """Write the cvt root metadata (detector_uids, number_of_simulated_events)."""
    detector_uids = Struct(
        {name: Scalar(int(uid)) for name, uid in zip(names, uids, strict=True)}
    )
    lh5.write(detector_uids, "detector_uids", str(path), wo_mode="append")
    lh5.write(Scalar(1000), "number_of_simulated_events", str(path), wo_mode="append")


def _trigger_table() -> Table:
    return Table(col_dict={"timestamp": Array(np.zeros(4, dtype=np.float32))})


def _geds_table() -> Table:
    """A geds sub-table with the single-template PSD fields the plot reads."""
    single_temp = Table(
        col_dict={
            "has_aoe": VectorOfVectors(data=[[True], [True], [True, True], [True]]),
            "is_single_site": VectorOfVectors(
                data=[[True], [False], [True, True], [True]]
            ),
        }
    )
    psd = Table(
        col_dict={
            "is_good": VectorOfVectors(data=[[True], [True], [True, True], [True]]),
            "single_temp": single_temp,
        }
    )
    return Table(
        col_dict={
            "energy": VectorOfVectors(
                data=[[500.0], [1000.0], [200.0, 300.0], [2000.0]]
            ),
            "is_good_channel": VectorOfVectors(
                data=[[True], [True], [True, True], [True]]
            ),
            "multiplicity": Array(np.array([1, 1, 2, 1], dtype=np.int32)),
            "psd": psd,
        }
    )


def _spms_table() -> Table:
    """An spms sub-table; ``time`` is doubly-jagged (per channel, per pe)."""
    return Table(
        col_dict={
            "multiplicity": Array(np.array([0, 2, 5, 1], dtype=np.int32)),
            "energy_sum": Array(np.array([0.0, 10.0, 30.0, 3.0], dtype=np.float32)),
            "time": VectorOfVectors(
                data=[[[]], [[100.0, 200.0], [150.0]], [[50.0]], [[300.0]]]
            ),
        }
    )


def _make_cvt_full(path: Path) -> None:
    """Full cvt file: both geds and spms tables, with coincident.spms."""
    coincident = Table(
        col_dict={
            "geds": Array(np.array([True, True, False, True])),
            "spms": Array(np.array([False, True, False, False])),
        }
    )
    evt = Table(
        col_dict={
            "trigger": _trigger_table(),
            "geds": _geds_table(),
            "spms": _spms_table(),
            "coincident": coincident,
        }
    )
    lh5.write(evt, "evt", str(path), wo_mode="write_safe")
    _write_cvt_root(path, ["V01", "B02"], [1, 2])


def _make_cvt_skip_opt(path: Path) -> None:
    """Skip-opt cvt file: no spms table, no coincident.spms."""
    coincident = Table(col_dict={"geds": Array(np.array([True, True, False, True]))})
    evt = Table(
        col_dict={
            "trigger": _trigger_table(),
            "geds": _geds_table(),
            "coincident": coincident,
        }
    )
    lh5.write(evt, "evt", str(path), wo_mode="write_safe")
    _write_cvt_root(path, ["V01", "B02"], [1, 2])


def _make_cvt_skip_hit(path: Path) -> None:
    """Skip-hit cvt file: no geds table, no coincident.geds."""
    coincident = Table(col_dict={"spms": Array(np.array([False, True, False, False]))})
    evt = Table(
        col_dict={
            "trigger": _trigger_table(),
            "spms": _spms_table(),
            "coincident": coincident,
        }
    )
    lh5.write(evt, "evt", str(path), wo_mode="write_safe")
    # skip_hit produces an empty detector_uids mapping at the evt tier
    _write_cvt_root(path, [], [])


def _run_plot(cvt_file: Path, out_pdf: Path, simid: str = "test_simid") -> None:
    """Execute the plot script with a mock ``snakemake`` global."""
    snakemake = SimpleNamespace(
        input=InputFiles(toclone=[str(cvt_file)]),
        output=[str(out_pdf)],
        wildcards=SimpleNamespace(simid=simid),
        config=AttrsDict({"nersc": {"dvs_ro": False}}),
    )
    try:
        runpy.run_path(
            str(_PLOT_SCRIPT),
            init_globals={"snakemake": snakemake},
            run_name="__main__",
        )
    finally:
        plt.close("all")


def test_plot_cvt_observables_full(tmp_path):
    """Both geds and spms present: the script builds all panels."""
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_full(cvt_file)
    out_pdf = tmp_path / "out.pdf"

    _run_plot(cvt_file, out_pdf)

    assert out_pdf.exists()
    assert out_pdf.stat().st_size > 0


def test_plot_cvt_observables_skip_opt(tmp_path):
    """No spms table (skip_opt): the script must skip the LAr/spms panels.

    Regression test: the plot used to access ``evt.coincident.spms`` and the
    ``spms`` sub-tables unconditionally, crashing with ``AttributeError`` when
    the evt tier was built with ``skip_opt``.
    """
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_skip_opt(cvt_file)
    out_pdf = tmp_path / "out.pdf"

    _run_plot(cvt_file, out_pdf)

    assert out_pdf.exists()
    assert out_pdf.stat().st_size > 0


def test_plot_cvt_observables_skip_hit(tmp_path):
    """No geds table (skip_hit): the script must skip the HPGe/energy panels."""
    cvt_file = tmp_path / "cvt.lh5"
    _make_cvt_skip_hit(cvt_file)
    out_pdf = tmp_path / "out.pdf"

    _run_plot(cvt_file, out_pdf)

    assert out_pdf.exists()
    assert out_pdf.stat().st_size > 0


@pytest.mark.needs_remage
def test_plot_cvt_observables_with_real_cvt(tmp_path, legend_cvt_path):
    """Run the plot on a real cvt file produced by the full test pipeline."""
    out_pdf = tmp_path / "out.pdf"

    _run_plot(legend_cvt_path, out_pdf, simid="birds_nest_K40")

    assert out_pdf.exists()
    assert out_pdf.stat().st_size > 0
