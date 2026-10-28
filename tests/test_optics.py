from __future__ import annotations

import logging
from pathlib import Path

import awkward as ak
import lh5
import numpy as np
import pytest
import reboost
import reboost.spms
from dbetto import AttrsDict
from lgdo import Array, Table, VectorOfVectors
from reboost.optmap.convolve import NumdetStats, OptmapForConvolve

from legendsimflow import optics

RUNID = "l200-p03-r001-phy"
OPTMAP = "remage/l200cfg01-optmap-dummy.lh5"

# the opt tier settings read by the optics functions
SETTINGS = AttrsDict(
    {
        "max_pes_per_hit_per_sipm": 0,
        "max_pes_per_hit_combined": 0,
        "store_expected_pes": False,
        "photoelectron_resolution_sigma": 0.0,
        "time_resolution_in_ns": 0.0,
    }
)


def _sipm_table(dt: list[list[float]], t0: float = 1000.0) -> Table:
    """A SiPM step table: photon arrival times `dt` after the hit time `t0`, in ns."""
    return Table(
        {
            "evtid": Array(np.arange(len(dt), dtype=np.int64)),
            "t0": Array(np.full(len(dt), t0), attrs={"units": "ns"}),
            "time": VectorOfVectors(ak.Array(dt) + t0, attrs={"units": "ns"}),
        }
    )


def _argon_table(n_hits: int, n_steps: int, seed: int = 42) -> Table:
    """An argon step table with `n_hits` hits of `n_steps` steps each.

    Steps are in the volume covered by the dummy optical map, in metres, with
    energy depositions of electrons in keV and times in ns.
    """
    rng = np.random.default_rng(seed)
    n = n_hits * n_steps
    t0 = np.arange(n_hits) * 1e6
    time = ak.unflatten(np.repeat(t0, n_steps) + rng.uniform(0, 10, n), n_steps)
    return Table(
        {
            "evtid": Array(np.arange(n_hits, dtype=np.int64)),
            "t0": Array(t0, attrs={"units": "ns"}),
            "time": VectorOfVectors(time),
            "edep": VectorOfVectors(ak.unflatten(rng.uniform(100, 1000, n), n_steps)),
            "particle": VectorOfVectors(
                ak.unflatten(np.full(n, 11, np.int32), n_steps)
            ),
            "xloc": VectorOfVectors(ak.unflatten(rng.uniform(-0.3, 0.3, n), n_steps)),
            "yloc": VectorOfVectors(ak.unflatten(rng.uniform(-0.3, 0.3, n), n_steps)),
            "zloc": VectorOfVectors(ak.unflatten(rng.uniform(0.5, 1.0, n), n_steps)),
        }
    )


def test_resolve_map_scaling_scalar():
    assert optics.resolve_map_scaling(0.3, "S002") == pytest.approx(0.3)
    assert optics.resolve_map_scaling(0.3, "all") == pytest.approx(0.3)
    assert optics.resolve_map_scaling(1, "S002") == pytest.approx(1.0)


def test_resolve_map_scaling_per_channel():
    setting = {"S002": 0.42, "default": 0.3}
    assert optics.resolve_map_scaling(setting, "S002") == pytest.approx(0.42)
    assert optics.resolve_map_scaling(setting, "S007") == pytest.approx(0.3)


def test_resolve_map_scaling_without_default_raises():
    with pytest.raises(KeyError, match="no optmap_scaling_factor entry for S007"):
        optics.resolve_map_scaling({"S002": 0.42}, "S007")


def test_get_usability(caplog):
    usability_map = {RUNID: {"S002": "ac"}}
    with caplog.at_level(logging.WARNING, logger="legendsimflow.optics"):
        assert optics.get_usability(usability_map, RUNID, "S002") == "ac"
        assert caplog.records == []
        assert optics.get_usability(usability_map, RUNID, "S003") == "on"
    assert [r.levelno for r in caplog.records] == [logging.WARNING]
    assert "S003" in caplog.records[0].message
    assert RUNID in caplog.records[0].message


def test_tracked_photons_to_pes():
    chunk = _sipm_table([[3.0, 1.0, 2.0], [5.0], []], t0=0).view_as("ak")

    pe_times, is_saturated = optics.tracked_photons_to_pes(chunk, 2)
    assert pe_times.to_list() == [[1.0, 2.0], [5.0], []]
    assert is_saturated.dtype == np.bool_
    assert is_saturated.tolist() == [True, False, False]

    for cap in (0, -1):
        pe_times, is_saturated = optics.tracked_photons_to_pes(chunk, cap)
        assert pe_times.to_list() == [[1.0, 2.0, 3.0], [5.0], []]
        assert is_saturated.dtype == np.bool_
        assert is_saturated.tolist() == [False, False, False]


def _constant_light(xloc, yloc, zloc, *, prob):  # noqa: ARG001
    return np.full_like(xloc, prob)


def test_load_light_model():
    model = optics.load_light_model(
        {"function": f"{__name__}._constant_light", "kwargs": {"prob": 0.2}}
    )
    assert model(*[np.array([100.0])] * 3) == [0.2]


def _geometry_light(xloc, yloc, zloc, *, geometry):  # noqa: ARG001
    return np.full_like(xloc, geometry)


def test_load_light_model_geometry():
    spec = {"function": f"{__name__}._geometry_light", "kwargs": {}}
    assert optics.load_light_model(spec, 0.1)(*[np.array([1.0])] * 3) == [0.1]
    with pytest.raises(ValueError, match="needs the GDML geometry"):
        optics.load_light_model(spec)
    # models without the argument do not get it
    spec = {"function": f"{__name__}._constant_light", "kwargs": {"prob": 0.2}}
    assert "geometry" not in optics.load_light_model(spec, 0.1).keywords


def test_model_photoelectrons():
    xyz = ak.Array([[0.0] * 3, [], [0.0]])
    photons = ak.Array([[10_000] * 3, [], [1]])

    n_pes, expected = optics.model_photoelectrons(
        lambda x, *_: np.full_like(x, 0.1),
        xyz,
        xyz,
        xyz,
        photons,
        max_pes_per_hit=1500,
    )
    assert ak.num(n_pes).to_list() == [3, 0, 1]
    assert expected == pytest.approx([3000, 0, 0.1])
    # about 1000 photoelectrons per step: the first two steps fill the hit up
    # to the cap
    assert ak.sum(n_pes[:1], axis=1).to_list() == [1500]
    assert 0 < n_pes[0, 1] < 1500
    assert n_pes[0, 2] == 0

    n_pes, expected = optics.model_photoelectrons(
        lambda x, *_: np.full_like(x, 0.1),
        xyz,
        xyz,
        xyz,
        photons,
        map_scaling=0,
    )
    assert ak.sum(n_pes) == 0
    assert expected == pytest.approx([3000, 0, 0.1])

    for bad in (
        lambda x, *_: x - 1,
        lambda x, *_: x + 2,
        lambda x, *_: np.full_like(x, np.nan),
        lambda *_: np.ones(2),
    ):
        with pytest.raises(ValueError, match="detection probability"):
            optics.model_photoelectrons(bad, xyz, xyz, xyz, photons)


def test_apply_photoelectron_resolution():
    pe_times = ak.Array([[0.0, 0.6, 1.1, 1.4, 2.3], [], [7.0]])

    # no resolution effects: times untouched, one photoelectron each
    times, amps = optics.apply_photoelectron_resolution(pe_times, 0.0, 0.0)
    assert times.to_list() == pe_times.to_list()
    assert amps.to_list() == [[1.0] * 5, [], [1.0]]

    # photoelectrons within 1 ns merge at the time of the first
    times, amps = optics.apply_photoelectron_resolution(pe_times, 0.0, 1.0)
    assert times.to_list() == [[0.0, 1.1, 2.3], [], [7.0]]
    assert amps.to_list() == [[2.0, 2.0, 1.0], [], [1.0]]

    # amplitude smearing: mean one, spread sigma
    pe_times = ak.Array([np.arange(10000, dtype=float) * 100])
    times, amps = optics.apply_photoelectron_resolution(pe_times, 0.1, 0.0)
    assert times.to_list() == pe_times.to_list()
    flat = ak.to_numpy(ak.flatten(amps))
    assert flat.mean() == pytest.approx(1.0, abs=0.01)
    assert flat.std() == pytest.approx(0.1, abs=0.01)


def test_build_opt_table():
    stp_chunk = _sipm_table([[1.0, 2.0], [], [5.0]], t0=1000.0)
    pe_times = ak.Array([[1001.0, 1002.0], [], [1005.0]])
    pe_amps = ak.Array([[1.0, 0.9], [], [1.1]])
    is_saturated = np.array([True, False, False])

    table = optics.build_opt_table(
        stp_chunk, pe_times, pe_amps, is_saturated, None, RUNID, "ac"
    )
    assert len(table) == 3
    assert set(table.keys()) == {
        "evtid", "t0", "dt", "energy", "is_saturated", "period", "run", "usability",
    }  # fmt: skip
    assert table.evtid.nda.tolist() == [0, 1, 2]
    assert table.t0.nda.tolist() == [1000.0] * 3
    assert table.dt.attrs["units"] == "ns"
    assert table.dt.flattened_data.nda.dtype == np.float32
    assert table.dt.view_as("ak").to_list() == [[1.0, 2.0], [], [5.0]]
    assert table.energy.flattened_data.nda.dtype == np.float32
    assert ak.flatten(table.energy.view_as("ak")).to_list() == pytest.approx(
        [1.0, 0.9, 1.1]
    )
    assert table.is_saturated.nda.dtype == np.bool_
    assert table.is_saturated.nda.tolist() == [True, False, False]
    assert table.period.nda.tolist() == [3] * 3
    assert table.run.nda.tolist() == [1] * 3
    assert table.usability.nda.tolist() == [1] * 3

    table = optics.build_opt_table(
        stp_chunk, pe_times, pe_amps, is_saturated, [0.5, 0.0, 2.0], RUNID, "on"
    )
    assert table.expected_pes.nda.dtype == np.float32
    assert table.expected_pes.nda.tolist() == [0.5, 0.0, 2.0]
    assert table.usability.nda.tolist() == [0] * 3


def test_build_opt_table_empty():
    stp_chunk = optics.empty_photon_table()
    pe_times, is_saturated = optics.tracked_photons_to_pes(stp_chunk.view_as("ak"), 2)
    pe_times, pe_amps = optics.apply_photoelectron_resolution(pe_times, 0.1, 10.0)
    table = optics.build_opt_table(
        stp_chunk, pe_times, pe_amps, is_saturated, None, RUNID, "on"
    )
    assert len(table) == 0
    assert "dt" in table
    assert table.dt.attrs["units"] == "ns"
    assert "period" in table


def test_empty_photon_table():
    table = optics.empty_photon_table()
    assert len(table) == 0
    assert set(table.keys()) == {"evtid", "t0", "time"}
    assert table.t0.attrs["units"] == "ns"
    assert table.time.attrs["units"] == "ns"
    assert table.time.view_as("ak").to_list() == []


def test_log_optmap_coverage(caplog):
    stats = NumdetStats(
        energy_looped=100.0, energy_oob=2.0, energy_no_stats=1.0, energy_zero=0.5
    )
    with caplog.at_level(logging.DEBUG, logger="legendsimflow.optics"):
        optics.log_optmap_coverage(stats, "S002", "liquid_argon")
    assert [r.levelno for r in caplog.records] == [logging.WARNING, logging.DEBUG]
    assert "3.00% of the energy deposited in liquid_argon is outside the map" in (
        caplog.records[0].message
    )
    assert "0.50%" in caplog.records[1].message
    assert "zero detection probability" in caplog.records[1].message

    # nothing processed: no division by zero, nothing to warn about
    caplog.clear()
    with caplog.at_level(logging.DEBUG, logger="legendsimflow.optics"):
        optics.log_optmap_coverage(NumdetStats(), "all", "liquid_argon")
    assert [r.levelno for r in caplog.records] == [logging.DEBUG, logging.DEBUG]
    assert all("0.00%" in r.message for r in caplog.records)


def test_write_tracked_photoelectrons(tmp_path):
    out_file = tmp_path / "opt.lh5"
    chunks = [
        _sipm_table([[3.0, 1.0, 2.0], []], t0=1000.0),
        _sipm_table([[5.0]], t0=2000.0),
    ]
    settings = AttrsDict(SETTINGS | {"max_pes_per_hit_per_sipm": 2})

    optics.write_tracked_photoelectrons(
        chunks, out_file, "S002", 1002, RUNID, "ac", settings=settings
    )

    out = lh5.read_as("hit/S002", out_file, library="ak")
    assert out.evtid.to_list() == [0, 1, 0]
    assert out.t0.to_list() == [1000.0, 1000.0, 2000.0]
    assert out.dt.to_list() == [[1.0, 2.0], [], [5.0]]
    assert out.energy.to_list() == [[1.0, 1.0], [], [1.0]]
    assert out.is_saturated.to_list() == [True, False, False]
    assert "expected_pes" not in out.fields
    assert set(out.usability.to_list()) == {1}
    assert reboost.get_remage_detector_uids(out_file, lh5_table="hit") == {1002: "S002"}

    # a SiPM without photons gets an empty table with the same layout
    optics.write_tracked_photoelectrons(
        [optics.empty_photon_table()], out_file, "S003", 1003, RUNID, "on",
        settings=settings,
    )  # fmt: skip
    assert lh5.read_n_rows("hit/S003", out_file) == 0
    assert set(lh5.ls(out_file, "hit/S003/")) == {f"hit/S003/{f}" for f in out.fields}
    assert reboost.get_remage_detector_uids(out_file, lh5_table="hit") == {
        1002: "S002",
        1003: "S003",
    }


def test_optmap_to_pes(legend_testdata):
    optmap = reboost.spms.load_optmap(legend_testdata[OPTMAP], "all")
    chunk = _argon_table(n_hits=50, n_steps=5).view_as("ak")

    pe_times, is_saturated, expected_pes, stats = optics.optmap_to_pes(
        chunk, optmap, "all", material="lar", map_scaling=1.0, max_pes_per_hit=0,
        store_expected_pes=True,
    )  # fmt: skip
    assert len(pe_times) == 50
    assert is_saturated.dtype == np.bool_
    assert not np.any(is_saturated)
    assert len(expected_pes) == 50
    assert np.all(np.asarray(expected_pes) > 0)
    assert ak.sum(ak.num(pe_times)) > 0
    # sorted and after the first step of the hit
    assert ak.all(pe_times[:, 1:] >= pe_times[:, :-1])
    assert ak.all(pe_times >= chunk.t0)
    assert isinstance(stats, NumdetStats)
    assert stats.steps_looped == 250
    assert stats.steps_oob == 0

    # a smaller scaling factor gives fewer photoelectrons, in proportion. the
    # expectation is at unit efficiency, so it does not change
    pe_times_half, _, expected_half, _ = optics.optmap_to_pes(
        chunk, optmap, "all", material="lar", map_scaling=0.5, max_pes_per_hit=0,
        store_expected_pes=True,
    )  # fmt: skip
    n_pes = ak.sum(ak.num(pe_times))
    assert n_pes > 10000
    assert ak.sum(ak.num(pe_times_half)) / n_pes == pytest.approx(0.5, rel=0.05)
    assert np.mean(expected_half) == pytest.approx(np.mean(expected_pes), rel=0.05)

    # the cap flags saturated hits and nothing is returned for the expectation
    pe_times, is_saturated, expected_pes, _ = optics.optmap_to_pes(
        chunk, optmap, "all", material="lar", map_scaling=1.0, max_pes_per_hit=1,
        store_expected_pes=False,
    )  # fmt: skip
    assert expected_pes is None
    assert ak.all(ak.num(pe_times) <= 1)
    assert is_saturated.tolist() == (ak.num(pe_times) >= 1).to_list()
    assert np.any(is_saturated)


def test_group_by_window():
    grouped = optics.group_by_window(
        ak.Array([[1.0], [2.0, 3.0], [], [4.0]]), [2, 0, 2, 0], 4
    )
    assert grouped.to_list() == [[[2.0, 3.0], [4.0]], [], [[1.0], []], []]
    assert optics.group_by_window(np.array([5, 6, 7]), [1, 1, 0], 2).to_list() == [
        [7],
        [5, 6],
    ]


def test_find_event_windows():
    # UIDs 1 and 2 are light sources, 3 is a germanium detector. the hits of
    # a window are ordered by time, not by table
    tcm = ak.Array(
        {
            "table_key": [[2, 1], [3], [2], [1, 1], [3, 2]],
            "row_in_table": [[5, 0], [0], [6], [2, 1], [1, 7]],
        }
    )
    windows, n_windows = optics.find_event_windows(tcm, {1: "pen", 2: "lar"})
    assert n_windows == 4
    assert windows.keys() == {"lar", "pen"}
    assert windows["lar"][0] == 5
    assert windows["lar"][1].tolist() == [0, 1, 3]
    assert windows["pen"][0] == 0
    assert windows["pen"][1].tolist() == [0, 2, 2]

    windows, n_windows = optics.find_event_windows(tcm[1:2], {1: "pen", 2: "lar"})
    assert (windows, n_windows) == ({}, 0)

    # a gap in the rows of a source
    tcm = ak.Array({"table_key": [[1], [1]], "row_in_table": [[0], [2]]})
    with pytest.raises(RuntimeError, match="pen are not one block"):
        optics.find_event_windows(tcm, {1: "pen"})


def test_pool_photoelectrons():
    pe_times = ak.Array([[3.0, 1.0], [2.0], [7.0], [0.5]])
    windows, pooled, is_saturated, expected = optics.pool_photoelectrons(
        pe_times,
        np.array([1.0, 2.0, 3.0, 4.0]),
        t0=[10.0, 5.0, 7.0, 0.0],
        evtid=[4, 4, 5, 6],
        window_idx=[0, 0, 1, 2],
        n_windows=3,
        max_pes_per_hit=2,
    )
    assert pooled.to_list() == [[1.0, 2.0], [7.0], [0.5]]
    assert is_saturated.tolist() == [True, False, False]
    assert expected.tolist() == [3.0, 3.0, 4.0]
    # the earliest hit of the window
    assert windows.t0.nda.tolist() == [5.0, 7.0, 0.0]
    assert windows.t0.attrs["units"] == "ns"
    assert windows.evtid.nda.tolist() == [4, 5, 6]

    _, pooled, is_saturated, expected = optics.pool_photoelectrons(
        pe_times, None, [0.0] * 4, [0] * 4, [0, 0, 1, 2], 3, 0
    )
    assert pooled.to_list() == [[1.0, 2.0, 3.0], [7.0], [0.5]]
    assert not np.any(is_saturated)
    assert expected is None


def _write_argon_file(path: Path, n_hits: int, seed: int = 42) -> Path:
    reboost.write_hit_table_chunk(
        _argon_table(n_hits, 5, seed=seed), "stp/liquid_argon", path, uid=1
    )
    return path


def _one_source(scaling: float | dict = 1.0) -> dict:
    return {
        "lar": AttrsDict(
            {
                "table": "liquid_argon",
                "material": "lar",
                "optmap_scaling_factor": scaling,
            }
        )
    }


def test_write_scintillation_photoelectrons(legend_testdata, tmp_path):
    stp_file = _write_argon_file(tmp_path / "stp.lh5", 30)
    out_file = tmp_path / "opt.lh5"
    optmap_file = legend_testdata[OPTMAP]
    sources = _one_source({"S002": 0.5, "default": 1.0})
    settings = AttrsDict(
        SETTINGS
        | {
            "max_pes_per_hit_combined": 1,
            "store_expected_pes": True,
            "time_resolution_in_ns": 10.0,
        }
    )

    # combined map, pre-loaded: table hit/spms, cap of the combined map. one
    # hit per window, read in two chunks
    optmap = reboost.spms.load_optmap(optmap_file, "all")
    stats = optics.write_scintillation_photoelectrons(
        stp_file, {"lar": (0, np.arange(30))}, 30, sources, {"lar": optmap},
        out_file, "all", 1, RUNID, "on", settings=settings, rows_per_chunk=20,
    )  # fmt: skip
    assert stats.keys() == {"lar"}
    assert stats["lar"].steps_looped == 150
    out = lh5.read_as("hit/spms", out_file, library="ak")
    assert len(out) == 30
    assert out.evtid.to_list() == list(range(30))
    assert ak.all(ak.num(out.dt) <= 1)
    assert ak.any(ak.num(out.dt) == 1)
    assert out.is_saturated.to_list() == (ak.num(out.dt) == 1).to_list()
    assert ak.all(out.expected_pes > 0)
    assert ak.all(out.dt >= 0)
    assert ak.all(out.energy == 1)

    # per-SiPM map loaded from file: table hit/<sipm>, cap of the single SiPM
    stp_file = _write_argon_file(tmp_path / "stp-large.lh5", 2000, seed=1)
    windows = {"lar": (0, np.arange(2000))}
    optics.write_scintillation_photoelectrons(
        stp_file, windows, 2000, sources, {"lar": optmap_file}, out_file, "S002",
        1002, RUNID, "ac", settings=settings, rows_per_chunk=500,
    )  # fmt: skip
    out = lh5.read_as("hit/S002", out_file, library="ak")
    assert len(out) == 2000
    assert not ak.any(out.is_saturated)
    assert set(out.usability.to_list()) == {1}
    assert reboost.get_remage_detector_uids(out_file, lh5_table="hit") == {
        1: "spms",
        1002: "S002",
    }

    # the per-channel scaling factor halves the number of photoelectrons. with
    # no amplitude smearing the summed amplitude counts them, also after the
    # time clustering
    optics.write_scintillation_photoelectrons(
        stp_file, windows, 2000, _one_source(), {"lar": optmap_file},
        tmp_path / "ref.lh5", "S002", 1002, RUNID, "ac", settings=settings,
        rows_per_chunk=500,
    )  # fmt: skip
    ref = lh5.read_as("hit/S002", tmp_path / "ref.lh5", library="ak")
    assert ak.sum(ref.energy) > 10000
    assert ak.sum(out.energy) / ak.sum(ref.energy) == pytest.approx(0.5, rel=0.1)

    # no window with scintillator hits: nothing to write
    stats = optics.write_scintillation_photoelectrons(
        stp_file, {}, 0, sources, {"lar": optmap_file}, tmp_path / "none.lh5",
        "S002", 1002, RUNID, "ac", settings=settings, rows_per_chunk=500,
    )  # fmt: skip
    assert stats == {}
    assert not (tmp_path / "none.lh5").exists()

    # a light source without hits: its map file is not read, here it does not
    # even exist
    pen = {
        "pen": AttrsDict(
            {"table": "pen", "material": "pen", "optmap_scaling_factor": 1.0}
        )
    }
    stats = optics.write_scintillation_photoelectrons(
        stp_file, windows, 2000, sources | pen,
        {"lar": optmap_file, "pen": tmp_path / "missing.lh5"}, tmp_path / "pen.lh5",
        "S002", 1002, RUNID, "ac", settings=settings, rows_per_chunk=500,
    )  # fmt: skip
    assert stats.keys() == {"lar"}


def test_write_scintillation_photoelectrons_pooled(legend_testdata, tmp_path):
    stp_file = _write_argon_file(tmp_path / "stp.lh5", 6)
    optmap = reboost.spms.load_optmap(legend_testdata[OPTMAP], "all")
    settings = AttrsDict(SETTINGS | {"store_expected_pes": True})

    def run(windows, n_windows, out_file):
        optics.write_scintillation_photoelectrons(
            stp_file, {"lar": (0, np.asarray(windows))}, n_windows, _one_source(),
            {"lar": optmap}, out_file, "all", 1, RUNID, "on", settings=settings,
            rows_per_chunk=1,
        )  # fmt: skip
        return lh5.read_as("hit/spms", out_file, library="ak")

    single = run(range(6), 6, tmp_path / "single.lh5")
    # pairs of hits in the same window: the windows are never split between
    # chunks, and their light adds up
    pooled = run([0, 0, 1, 1, 2, 2], 3, tmp_path / "pooled.lh5")
    assert len(pooled) == 3
    assert pooled.evtid.to_list() == [0, 2, 4]
    assert pooled.t0.to_list() == single.t0[::2].to_list()
    assert ak.all(pooled.dt[:, 1:] >= pooled.dt[:, :-1])
    # the number of emitted photons is sampled, hence the tolerance
    expected = ak.to_numpy(single.expected_pes)
    assert ak.to_numpy(pooled.expected_pes) == pytest.approx(
        expected[::2] + expected[1::2], rel=0.05
    )


def test_write_scintillation_photoelectrons_light_model(legend_testdata, tmp_path):
    stp_file = _write_argon_file(tmp_path / "stp.lh5", 4)
    reboost.write_hit_table_chunk(
        _argon_table(4, 5, seed=1), "stp/pen", stp_file, uid=2
    )
    sources = _one_source() | {
        "pen": AttrsDict(
            {"table": "pen", "material": "pen", "optmap_scaling_factor": 1.0}
        )
    }
    model = optics.load_light_model(
        {"function": f"{__name__}._constant_light", "kwargs": {"prob": 0.01}}
    )
    settings = AttrsDict(SETTINGS | {"store_expected_pes": True})

    # one argon and one PEN hit per window
    windows = {"lar": (0, np.arange(4)), "pen": (0, np.arange(4))}
    stats = optics.write_scintillation_photoelectrons(
        stp_file, windows, 4, sources,
        {"lar": legend_testdata[OPTMAP], "pen": model}, tmp_path / "opt.lh5",
        "all", 1, RUNID, "on", settings=settings, rows_per_chunk=4,
    )  # fmt: skip
    # no coverage statistics for a light model
    assert stats.keys() == {"lar"}
    out = lh5.read_as("hit/spms", tmp_path / "opt.lh5", library="ak")
    assert len(out) == 4
    assert ak.all(out.expected_pes > 0)
    assert ak.all(out.dt[:, 1:] >= out.dt[:, :-1])


def test_write_photoelectrons_profiler(legend_testdata, tmp_path):
    """The reboost profiler can time the blocks of the SiPM processing."""
    perf_block, print_perf, _ = reboost.make_profiler()
    optics.write_tracked_photoelectrons(
        [_sipm_table([[1.0]])], tmp_path / "opt.lh5", "S002", 1002, RUNID, "on",
        settings=SETTINGS, perf_block=perf_block,
    )  # fmt: skip
    optics.write_scintillation_photoelectrons(
        _write_argon_file(tmp_path / "stp.lh5", 3), {"lar": (0, np.arange(3))}, 3,
        _one_source(), {"lar": legend_testdata[OPTMAP]}, tmp_path / "opt.lh5",
        "all", 1, RUNID, "on", settings=SETTINGS, rows_per_chunk=2,
        perf_block=perf_block,
    )  # fmt: skip
    print_perf()
    assert {"hit/S002", "hit/spms"} <= set(lh5.ls(tmp_path / "opt.lh5", "hit/*"))


def test_optmap_object_is_kept(legend_testdata, tmp_path, monkeypatch):
    """A map given as an object is not reloaded from disk."""
    optmap = reboost.spms.load_optmap(legend_testdata[OPTMAP], "all")
    assert isinstance(optmap, OptmapForConvolve)
    monkeypatch.setattr(reboost.spms, "load_optmap", lambda *_: pytest.fail("reloaded"))
    optics.write_scintillation_photoelectrons(
        _write_argon_file(tmp_path / "stp.lh5", 5), {"lar": (0, np.arange(5))}, 5,
        _one_source(), {"lar": optmap}, tmp_path / "opt.lh5", "all", 1, RUNID, "on",
        settings=SETTINGS, rows_per_chunk=5,
    )  # fmt: skip
