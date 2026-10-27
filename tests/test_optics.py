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
        "light_source": "tracked_photons",
        "optmap_scaling_factor": 1.0,
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


def test_write_sipm_photoelectrons_tracked_photons(tmp_path):
    out_file = tmp_path / "opt.lh5"
    chunks = [
        _sipm_table([[3.0, 1.0, 2.0], []], t0=1000.0),
        _sipm_table([[5.0]], t0=2000.0),
    ]
    settings = AttrsDict(SETTINGS | {"max_pes_per_hit_per_sipm": 2})

    stats = optics.write_sipm_photoelectrons(
        chunks, out_file, "S002", 1002, RUNID, "ac", settings=settings
    )
    assert stats is None

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
    optics.write_sipm_photoelectrons(
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
        chunk, optmap, "all", map_scaling=1.0, max_pes_per_hit=0,
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
        chunk, optmap, "all", map_scaling=0.5, max_pes_per_hit=0,
        store_expected_pes=True,
    )  # fmt: skip
    n_pes = ak.sum(ak.num(pe_times))
    assert n_pes > 10000
    assert ak.sum(ak.num(pe_times_half)) / n_pes == pytest.approx(0.5, rel=0.05)
    assert np.mean(expected_half) == pytest.approx(np.mean(expected_pes), rel=0.05)

    # the cap flags saturated hits and nothing is returned for the expectation
    pe_times, is_saturated, expected_pes, _ = optics.optmap_to_pes(
        chunk, optmap, "all", map_scaling=1.0, max_pes_per_hit=1,
        store_expected_pes=False,
    )  # fmt: skip
    assert expected_pes is None
    assert ak.all(ak.num(pe_times) <= 1)
    assert is_saturated.tolist() == (ak.num(pe_times) >= 1).to_list()
    assert np.any(is_saturated)


def test_write_sipm_photoelectrons_optmap(legend_testdata, tmp_path):
    out_file = tmp_path / "opt.lh5"
    optmap_file = legend_testdata[OPTMAP]
    chunks = [_argon_table(20, 5, seed=1), _argon_table(10, 5, seed=2)]
    settings = AttrsDict(
        SETTINGS
        | {
            "light_source": "optmap",
            "optmap_scaling_factor": {"S002": 0.5, "default": 1.0},
            "max_pes_per_hit_combined": 1,
            "store_expected_pes": True,
            "time_resolution_in_ns": 10.0,
        }
    )

    # combined map, pre-loaded: table hit/spms, cap of the combined map
    optmap = reboost.spms.load_optmap(optmap_file, "all")
    stats = optics.write_sipm_photoelectrons(
        chunks, out_file, "all", 1, RUNID, "on", settings=settings, optmap=optmap
    )
    assert isinstance(stats, NumdetStats)
    assert stats.steps_looped == 150
    out = lh5.read_as("hit/spms", out_file, library="ak")
    assert len(out) == 30
    assert out.evtid.to_list() == [*range(20), *range(10)]
    assert ak.all(ak.num(out.dt) <= 1)
    assert ak.any(ak.num(out.dt) == 1)
    assert out.is_saturated.to_list() == (ak.num(out.dt) == 1).to_list()
    assert ak.all(out.expected_pes > 0)
    assert ak.all(out.dt >= 0)
    assert ak.all(out.energy == 1)

    # per-SiPM map loaded from file: table hit/<sipm>, cap of the single SiPM
    stats = optics.write_sipm_photoelectrons(
        [_argon_table(2000, 5, seed=1)], out_file, "S002", 1002, RUNID, "ac",
        settings=settings, optmap=optmap_file,
    )  # fmt: skip
    assert isinstance(stats, NumdetStats)
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
    settings.optmap_scaling_factor = 1.0
    optics.write_sipm_photoelectrons(
        [_argon_table(2000, 5, seed=1)], tmp_path / "ref.lh5", "S002", 1002, RUNID,
        "ac", settings=settings, optmap=optmap_file,
    )  # fmt: skip
    ref = lh5.read_as("hit/S002", tmp_path / "ref.lh5", library="ak")
    assert ak.sum(ref.energy) > 10000
    assert ak.sum(out.energy) / ak.sum(ref.energy) == pytest.approx(0.5, rel=0.1)


def test_write_sipm_photoelectrons_profiler(tmp_path):
    """The reboost profiler can time the blocks of the SiPM processing."""
    perf_block, print_perf, _ = reboost.make_profiler()
    optics.write_sipm_photoelectrons(
        [_sipm_table([[1.0]])], tmp_path / "opt.lh5", "S002", 1002, RUNID, "on",
        settings=SETTINGS, perf_block=perf_block,
    )  # fmt: skip
    print_perf()
    assert Path(tmp_path / "opt.lh5").exists()


def test_optmap_object_is_kept(legend_testdata, tmp_path, monkeypatch):
    """A map given as an object is not reloaded from disk."""
    optmap = reboost.spms.load_optmap(legend_testdata[OPTMAP], "all")
    assert isinstance(optmap, OptmapForConvolve)
    monkeypatch.setattr(reboost.spms, "load_optmap", lambda *_: pytest.fail("reloaded"))
    optics.write_sipm_photoelectrons(
        [_argon_table(5, 3)], tmp_path / "opt.lh5", "all", 1, RUNID, "on",
        settings=AttrsDict(SETTINGS | {"light_source": "optmap"}), optmap=optmap,
    )  # fmt: skip
