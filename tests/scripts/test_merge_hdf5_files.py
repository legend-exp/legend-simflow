from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pytest

SCRIPT = (
    Path(__file__).parents[2] / "workflow/src/legendsimflow/scripts/merge_hdf5_files.sh"
)

pytestmark = pytest.mark.skipif(
    shutil.which("h5copy") is None, reason="HDF5 command-line tools not found"
)


def _run(*args):
    return subprocess.run(
        ["bash", str(SCRIPT), *map(str, args)],
        capture_output=True,
        text=True,
        check=False,
    )


def test_merge(tmp_path):
    with h5py.File(tmp_path / "a.h5", "w") as f:
        f.attrs["root"] = "a"
        f.create_group("V01").attrs["datatype"] = "struct{x}"
        f["V01/x"] = np.arange(3)
    with h5py.File(tmp_path / "b.h5", "w") as f:
        f.create_group("V 02/sub")
        f["V 02/sub/x"] = np.arange(5)
        f["V 02"].attrs["datatype"] = "struct{sub}"
        f["V03"] = 1.5

    res = _run(tmp_path / "out.h5", tmp_path / "a.h5", tmp_path / "b.h5")
    assert res.returncode == 0, res.stderr

    with h5py.File(tmp_path / "out.h5", "r") as f:
        assert set(f.keys()) == {"V01", "V 02", "V03"}
        assert f.attrs["root"] == "a"
        assert f["V01"].attrs["datatype"] == "struct{x}"
        assert f["V 02"].attrs["datatype"] == "struct{sub}"
        np.testing.assert_array_equal(f["V01/x"][:], np.arange(3))
        np.testing.assert_array_equal(f["V 02/sub/x"][:], np.arange(5))
        assert f["V03"][()] == 1.5


def test_merge_no_inputs(tmp_path):
    res = _run(tmp_path / "out.h5")
    assert res.returncode == 0, res.stderr

    with h5py.File(tmp_path / "out.h5", "r") as f:
        assert len(f.keys()) == 0


def test_merge_duplicate_name(tmp_path):
    for name in ("a.h5", "b.h5"):
        with h5py.File(tmp_path / name, "w") as f:
            f["V01"] = 1

    res = _run(tmp_path / "out.h5", tmp_path / "a.h5", tmp_path / "b.h5")
    assert res.returncode != 0
    assert "'/V01'" in res.stderr


def test_usage():
    assert _run().returncode != 0
