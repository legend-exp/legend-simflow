from __future__ import annotations

import logging
import sys
import warnings
from pathlib import Path

import awkward as ak
import lh5
import numpy as np
import pyg4ometry
import pygeomtools
import pytest
import reboost
import yaml
from dbetto import AttrsDict
from lgdo import Array, Table, VectorOfVectors

from legendsimflow.exceptions import SimflowConfigError
from legendsimflow.scripts.tier import opt
from legendsimflow.tcm import merge_stp_n_opt_tcms

dummyprod = Path(__file__).parent.parent / "dummyprod"

# channels with non-trivial usabilities added to the dummyprod status file
_SIPM_ON = "S002"  # usability: on  → code 0
_SIPM_AC = "S003"  # usability: ac  → code 1
_SIPM_OFF = "S007"  # usability: off → code 2


@pytest.mark.needs_remage
def test_opt_script_cli(
    tmp_path,
    monkeypatch,
    legend_testdata,
    legend_stp_path,
    legend_gdml_path,
    legend_simstat_part_path,
    legend_detector_usabilities_path,
):
    optmap_path = Path(legend_testdata.get_path("remage/l200cfg01-optmap-dummy.lh5"))

    raw = yaml.safe_load((dummyprod / "simflow-config-l200.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    config_path = tmp_path / "simflow-config-l200.yaml"
    config_path.write_text(yaml.safe_dump(raw))

    common_args = [
        "--stp-file", str(legend_stp_path),
        "--optmap", f"lar={optmap_path}",
        "--geom-file", str(legend_gdml_path),
        "--usability-file", str(legend_detector_usabilities_path / "usability.yaml"),
        "--jobid", "0000",
        "--simflow-config", str(config_path),
    ]  # fmt: skip

    # --- run 1: sum optical map (all SiPMs together) ---
    opt_file = tmp_path / "opt.lh5"
    monkeypatch.setattr(sys, "argv", [
        "opt", *common_args,
        "--simstat-part-file", str(legend_simstat_part_path),
        "--opt-file", str(opt_file),
    ])  # fmt: skip
    opt.main()

    assert opt_file.exists(), "opt output file was not created"

    root_keys = lh5.ls(opt_file)
    assert "hit" in root_keys, f"'hit' group missing; got {root_keys}"
    assert "tcm" in root_keys, f"'tcm' group missing; got {root_keys}"
    assert "hit/spms" in lh5.ls(opt_file, "hit/"), (
        f"'hit/spms' missing; got {lh5.ls(opt_file, 'hit/')}"
    )

    expected_fields = {
        "dt", "energy", "evtid", "is_saturated", "period", "run", "t0", "usability",
    }  # fmt: skip
    spms_fields = {f.removeprefix("hit/spms/") for f in lh5.ls(opt_file, "hit/spms/")}
    assert not (expected_fields - spms_fields), (
        f"fields {expected_fields - spms_fields} missing from hit/spms"
    )
    # renamed to dt: readers of the old field must fail
    assert "time" not in spms_fields, "legacy field 'time' still in hit/spms"
    assert "expected_pes" not in spms_fields, "expected_pes stored without the setting"

    def _field(table: str, name: str, f: Path = opt_file) -> np.ndarray:
        return lh5.read_as(f"hit/{table}/{name}", f, library="np")

    assert np.all(_field("spms", "period") == 3)
    assert np.all(np.isin(_field("spms", "run"), [0, 1]))
    assert np.all(np.isin(_field("spms", "usability"), [0, 1, 2]))

    energy = _field("spms", "energy")
    assert np.any(np.isfinite(energy)), "energy contains no finite values"
    assert np.all(energy[np.isfinite(energy)] >= 0), "energy contains negative values"

    is_sat = _field("spms", "is_saturated")
    assert is_sat.dtype == np.bool_, f"unexpected is_saturated dtype {is_sat.dtype}"

    # relative to the hit t0: a few microseconds at most
    dt = lh5.read("hit/spms/dt", opt_file).flattened_data.nda
    assert len(dt) > 0, "no photoelectrons in the opt file"
    assert np.all(np.abs(dt) < 1e5), f"dt is not relative to t0: max {np.max(dt)} ns"

    # --- run 2: per-SiPM mode with tiny partition → check usability codes ---
    part_file = tmp_path / "partitions_small.yaml"
    part_file.write_text(yaml.safe_dump({"job_0000": {"l200-p03-r000-phy": [0, 99]}}))

    opt_per_sipm = tmp_path / "opt_per_sipm.lh5"
    monkeypatch.setattr(sys, "argv", [
        "opt", *common_args,
        "--simstat-part-file", str(part_file),
        "--opt-file", str(opt_per_sipm),
        "--optmap-per-sipm",
    ])  # fmt: skip
    settings = opt.get_tier_settings
    monkeypatch.setattr(
        opt,
        "get_tier_settings",
        lambda config, tier: AttrsDict(
            settings(config, tier) | {"store_expected_pes": True}
        ),
    )
    opt.main()

    hit_tables = {t.removeprefix("hit/") for t in lh5.ls(opt_per_sipm, "hit/")}
    for sipm in (_SIPM_ON, _SIPM_AC, _SIPM_OFF):
        assert sipm in hit_tables, f"expected hit/{sipm} not found; got {hit_tables}"
        fields = {
            f.removeprefix(f"hit/{sipm}/") for f in lh5.ls(opt_per_sipm, f"hit/{sipm}/")
        }
        assert "dt" in fields, f"'dt' missing from hit/{sipm}"
        assert "time" not in fields, f"legacy field 'time' still in hit/{sipm}"
        assert "expected_pes" in fields, f"'expected_pes' missing from hit/{sipm}"
        assert np.all(_field(sipm, "expected_pes", opt_per_sipm) >= 0)

    assert np.all(_field(_SIPM_ON, "usability", opt_per_sipm) == 0), (
        f"{_SIPM_ON} usability should be 0 (on)"
    )
    assert np.all(_field(_SIPM_AC, "usability", opt_per_sipm) == 1), (
        f"{_SIPM_AC} usability should be 1 (ac)"
    )
    assert np.all(_field(_SIPM_OFF, "usability", opt_per_sipm) == 2), (
        f"{_SIPM_OFF} usability should be 2 (off)"
    )


@pytest.mark.needs_remage
def test_opt_script_tracked_photons(
    tmp_path,
    monkeypatch,
    legend_tracked_photons_stp_path,
    legend_gdml_path,
    legend_detector_usabilities_path,
):
    stp_file = legend_tracked_photons_stp_path
    # the argon table is only needed with an optical map
    assert "stp/liquid_argon" not in lh5.ls(stp_file, "stp/")
    photons = lh5.read_as(f"stp/{_SIPM_ON}", stp_file, library="ak")
    evtid = ak.to_numpy(photons.evtid)
    t0 = ak.to_numpy(photons.t0)
    n_ph = ak.to_numpy(ak.num(photons.time))
    sens_tables = pygeomtools.detectors.get_all_senstables(
        pyg4ometry.gdml.Reader(str(legend_gdml_path)).getRegistry()
    )

    n_tcm = lh5.read_n_rows("tcm", stp_file)
    part_file = tmp_path / "partitions.yaml"
    part_file.write_text(
        yaml.safe_dump({"job_0000": {"l200-p03-r000-phy": [0, n_tcm - 1]}})
    )

    raw = yaml.safe_load((dummyprod / "simflow-config-l200.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    config_path = tmp_path / "simflow-config-l200.yaml"
    config_path.write_text(yaml.safe_dump(raw))

    opt_file = tmp_path / "opt.lh5"
    monkeypatch.setattr(sys, "argv", [
        "opt",
        "--stp-file", str(stp_file),
        "--geom-file", str(legend_gdml_path),
        "--simstat-part-file", str(part_file),
        "--usability-file", str(legend_detector_usabilities_path / "usability.yaml"),
        "--jobid", "0000",
        "--opt-file", str(opt_file),
        "--simflow-config", str(config_path),
    ])  # fmt: skip
    settings = opt.get_tier_settings
    monkeypatch.setattr(
        opt,
        "get_tier_settings",
        lambda config, tier: AttrsDict(
            settings(config, tier)
            | {
                "light_source": "tracked_photons",
                # ignored with tracked photons
                "optmap_scaling_factor": 0.3,
                "store_expected_pes": True,
                "max_pes_per_hit_per_sipm": 2,
            }
        ),
    )
    opt.main()

    # one row per row of the SiPM table, also the one without argon
    out = lh5.read_as(f"hit/{_SIPM_ON}", opt_file, library="ak")
    assert out.evtid.to_list() == evtid.tolist()
    assert np.allclose(ak.to_numpy(out.t0), t0)
    assert "expected_pes" not in out.fields
    # capped at 2 photoelectrons per row, the first ones in time
    n_kept = np.minimum(n_ph, 2)
    assert ak.num(out.dt).to_list() == n_kept.tolist()
    expected_dt = ak.Array([1000 * np.arange(n) for n in n_kept])
    assert ak.all(abs(ak.flatten(out.dt - expected_dt)) < 1e-2)
    assert out.is_saturated.to_list() == (n_ph >= 2).tolist()

    # a SiPM without photons gets an empty table
    assert lh5.read_n_rows(f"hit/{_SIPM_AC}", opt_file) == 0
    assert sens_tables[_SIPM_AC].uid in reboost.get_remage_detector_uids(
        opt_file, lh5_table="hit"
    )


def _write_config(tmp_path: Path) -> Path:
    raw = yaml.safe_load((dummyprod / "simflow-config-l200.yaml").read_text())
    raw["paths"]["metadata"] = str(dummyprod / "inputs")
    config_path = tmp_path / "simflow-config-l200.yaml"
    config_path.write_text(yaml.safe_dump(raw))
    return config_path


def _override_settings(monkeypatch, **overrides) -> None:
    settings = opt.get_tier_settings
    monkeypatch.setattr(
        opt,
        "get_tier_settings",
        lambda config, tier: AttrsDict(settings(config, tier) | overrides),
    )


def _argv(tmp_path: Path, usability_dir: Path, **paths: str) -> list[str]:
    """Arguments of the opt script; `paths` overrides the input files."""
    files = {
        "stp-file": "stp.lh5",
        "optmap": "lar=optmap.lh5",
        "geom-file": "geom.gdml",
        "simstat-part-file": "partitions.yaml",
    } | {k.replace("_", "-"): v for k, v in paths.items()}
    return [
        "opt",
        *[a for k, v in files.items() for a in (f"--{k}", str(v))],
        "--usability-file", str(usability_dir / "usability.yaml"),
        "--jobid", "0000",
        "--opt-file", str(tmp_path / "opt.lh5"),
        "--simflow-config", str(_write_config(tmp_path)),
    ]  # fmt: skip


def test_opt_script_invalid_light_source(
    tmp_path, monkeypatch, legend_detector_usabilities_path
):
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, legend_detector_usabilities_path))
    _override_settings(monkeypatch, light_source="candles")
    with pytest.raises(SimflowConfigError, match="light_source must be"):
        opt.main()


def test_opt_script_optmap_requires_map_file(
    tmp_path, monkeypatch, legend_detector_usabilities_path
):
    argv = _argv(tmp_path, legend_detector_usabilities_path)
    i = argv.index("--optmap")
    del argv[i : i + 2]
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SimflowConfigError, match="one optical map per light source"):
        opt.main()


@pytest.mark.needs_remage
def test_opt_script_unknown_scintillator_volume(
    tmp_path, monkeypatch, legend_gdml_path, legend_detector_usabilities_path
):
    argv = _argv(
        tmp_path, legend_detector_usabilities_path, geom_file=str(legend_gdml_path)
    )
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(
        opt,
        "get_opt_scintillators",
        lambda _: AttrsDict(
            {"lar": {"table": "not_a_volume", "material": "lar", "model": None}}
        ),
    )
    with pytest.raises(SimflowConfigError, match=r"not_a_volume.*liquid_argon"):
        opt.main()


@pytest.mark.needs_remage
def test_opt_script_missing_scintillator_table(
    tmp_path,
    monkeypatch,
    caplog,
    legend_testdata,
    legend_tracked_photons_stp_path,
    legend_gdml_path,
    legend_detector_usabilities_path,
):
    """With an optical map, no argon table gives a warning and no opt file.

    Tracked photons do not need the table, see
    :func:`test_opt_script_tracked_photons`.
    """
    stp_file = legend_tracked_photons_stp_path
    assert "stp/liquid_argon" not in lh5.ls(stp_file, "stp/")
    n_tcm = lh5.read_n_rows("tcm", stp_file)
    part_file = tmp_path / "partitions.yaml"
    part_file.write_text(
        yaml.safe_dump({"job_0000": {"l200-p03-r000-phy": [0, n_tcm - 1]}})
    )
    monkeypatch.setattr(sys, "argv", _argv(
        tmp_path, legend_detector_usabilities_path,
        stp_file=str(stp_file),
        optmap=f"lar={legend_testdata.get_path('remage/l200cfg01-optmap-dummy.lh5')}",
        geom_file=str(legend_gdml_path),
        simstat_part_file=str(part_file),
    ))  # fmt: skip

    # nothing is written, so the TCM cannot be built at the end
    with caplog.at_level(logging.WARNING), pytest.raises(FileNotFoundError):
        opt.main()
    assert not (tmp_path / "opt.lh5").exists()
    assert any(
        "table liquid_argon of light source lar not found" in r.message
        and r.levelno == logging.WARNING
        for r in caplog.records
    )


@pytest.fixture
def two_sources_stp(tmp_path):
    """A box of argon holding two PEN plates and a SiPM, with hand-made hits.

    The PEN plates share the ``pen`` table. Event 0 has argon and PEN light in
    the same window, event 1 argon only, event 2 two PEN hits (written out of
    time order), event 3 argon and PEN 1 ms apart (two windows).
    """
    reg = pyg4ometry.geant4.Registry()
    world = pyg4ometry.geant4.LogicalVolume(
        pyg4ometry.geant4.solid.Box("world", 2, 2, 2, reg, lunit="m"),
        "G4_Galactic",
        "world",
        reg,
    )
    reg.setWorld(world)
    lar = pyg4ometry.geant4.LogicalVolume(
        pyg4ometry.geant4.solid.Box("lar", 1.5, 1.5, 1.5, reg, lunit="m"),
        "G4_lAr",
        "lar",
        reg,
    )
    pv = pyg4ometry.geant4.PhysicalVolume(
        [0, 0, 0], [0, 0, 0], lar, "liquid_argon", world, reg
    )
    pv.set_pygeom_active_detector(pygeomtools.RemageDetectorInfo("scintillator", 2))
    for i, name in enumerate(["pen_plate_1", "pen_plate_2", _SIPM_ON]):
        is_sipm = name == _SIPM_ON
        lv = pyg4ometry.geant4.LogicalVolume(
            pyg4ometry.geant4.solid.Box(name, 0.1, 0.1, 0.01, reg, lunit="m"),
            "G4_POLYETHYLENE",
            name,
            reg,
        )
        pv = pyg4ometry.geant4.PhysicalVolume(
            [0, 0, 0], [0, 0, 0.1 * (i + 1), "m"], lv, name, lar, reg
        )
        pv.set_pygeom_active_detector(
            pygeomtools.RemageDetectorInfo("optical", 3)
            if is_sipm
            else pygeomtools.RemageDetectorInfo(
                "scintillator", 1, allow_uid_reuse=True, ntuple_name="pen"
            )
        )
    gdml_file = tmp_path / "geom.gdml"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "optical detector .* no optical surface")
        pygeomtools.write_pygeom(reg, gdml_file)

    stp_file = tmp_path / "stp.lh5"
    hits = {
        # volume: (uid, edep per step in keV, [(evtid, t0 in ns), ...])
        "liquid_argon": (2, 10, [(0, 0), (1, 0), (3, 0)]),
        "pen": (1, 100, [(0, 100), (2, 60), (2, 50), (3, 1e6)]),
    }
    for volume, (uid, edep, rows) in hits.items():
        n = len(rows)

        def vov(value, units, n=n):
            return VectorOfVectors(ak.Array([[value] * 3] * n), attrs={"units": units})

        table = Table(
            {
                "evtid": Array(np.array([r[0] for r in rows])),
                "t0": Array(
                    np.array([r[1] for r in rows], dtype=float), attrs={"units": "ns"}
                ),
                "edep": vov(float(edep), "keV"),
                "particle": VectorOfVectors(ak.Array([[11] * 3] * n)),
                "xloc": vov(0.1, "m"),
                "yloc": vov(0.1, "m"),
                "zloc": vov(0.1, "m"),
                "time": VectorOfVectors(
                    ak.Array([[r[1], r[1] + 1, r[1] + 2] for r in rows]),
                    attrs={"units": "ns"},
                ),
            }
        )
        reboost.write_hit_table_chunk(table, f"stp/{volume}", stp_file, uid=uid)
    reboost.build_remage_tcm(stp_file, stp_file, lh5_group="stp")

    part_file = tmp_path / "partitions.yaml"
    part_file.write_text(yaml.safe_dump({"job_0000": {"l200-p03-r000-phy": [0, 4]}}))
    usability_file = tmp_path / "usability.yaml"
    usability_file.write_text(yaml.safe_dump({"l200-p03-r000-phy": {}}))

    return stp_file, gdml_file, part_file, usability_file


def _run_opt(monkeypatch, config_path, stp, optmaps, opt_file, per_sipm=False):
    stp_file, gdml_file, part_file, usability_file = stp
    monkeypatch.setattr(sys, "argv", [
        "opt",
        "--stp-file", str(stp_file),
        *[a for name, f in optmaps.items() for a in ("--optmap", f"{name}={f}")],
        "--geom-file", str(gdml_file),
        "--simstat-part-file", str(part_file),
        "--usability-file", str(usability_file),
        "--jobid", "0000",
        "--opt-file", str(opt_file),
        "--simflow-config", str(config_path),
        *(["--optmap-per-sipm"] if per_sipm else []),
    ])  # fmt: skip
    opt.main()
    table = _SIPM_ON if per_sipm else "spms"
    return lh5.read_as(f"hit/{table}", opt_file, library="ak")


def test_opt_script_two_light_sources(
    tmp_path, monkeypatch, legend_testdata, l200_config_factory, two_sources_stp
):
    optmap = legend_testdata.get_path("remage/l200cfg01-optmap-dummy.lh5")
    settings = {
        "store_expected_pes": True,
        "max_pes_per_hit_combined": 10_000,
        # two rows per chunk: windows are split across chunks
        "buffer_len": 2,
    }

    # without a scintillators block, only argon gives light
    config_path = l200_config_factory(tmp_path / "lar", {"opt": settings})
    out = _run_opt(
        monkeypatch, config_path, two_sources_stp, {"lar": optmap}, tmp_path / "lar.lh5"
    )
    assert out.evtid.to_list() == [0, 1, 3]

    config_path = l200_config_factory(
        tmp_path / "both",
        {
            "opt": settings
            | {
                "scintillators": {
                    "lar": {"table": "liquid_argon", "material": "lar"},
                    "pen": {
                        "table": "pen",
                        "material": "pen",
                        "optmap_scaling_factor": 0.5,
                    },
                }
            }
        },
    )
    opt_file = tmp_path / "both.lh5"
    out = _run_opt(
        monkeypatch,
        config_path,
        two_sources_stp,
        {"lar": optmap, "pen": optmap},
        opt_file,
    )

    # one row per event window with any scintillator hit
    assert out.evtid.to_list() == [0, 1, 2, 3, 3]
    assert ak.to_numpy(out.t0).tolist() == [0, 0, 50, 0, 1e6]
    assert ak.all(out.expected_pes > 0)
    # the smallest UID of the light sources (PEN), although lar comes first
    assert reboost.get_remage_detector_uids(opt_file, lh5_table="hit") == {1: "spms"}
    # the windows with PEN light only
    assert ak.all(ak.sum(out.energy[[2, 4]], axis=-1) > 0)
    # pooled photoelectrons are sorted in time and come after the window t0
    assert ak.all(out.dt[:, 1:] >= out.dt[:, :-1])
    assert ak.all(ak.flatten(out.dt) >= 0)

    # the evt tier puts each window back in place of all its scintillator hits
    stp_file = two_sources_stp[0]
    merged = merge_stp_n_opt_tcms(
        lh5.read_as("tcm", stp_file, library="ak"),
        lh5.read_as("tcm", opt_file, library="ak"),
        scintillator_uids=[1, 2],
    )
    assert ak.num(merged.table_key).to_list() == [1] * 5
    assert ak.flatten(merged.row_in_table).to_list() == [0, 1, 2, 3, 4]


def test_opt_script_two_light_sources_saturation(
    tmp_path, monkeypatch, legend_testdata, l200_config_factory, two_sources_stp
):
    optmap = legend_testdata.get_path("remage/l200cfg01-optmap-dummy.lh5")
    sources = {
        "lar": {"table": "liquid_argon", "material": "lar"},
        "pen": {"table": "pen", "material": "pen"},
    }
    config_path = l200_config_factory(
        tmp_path,
        {
            "opt": {
                # windows have about 20 (argon only) to 100 photoelectrons
                "max_pes_per_hit_combined": 40,
                # no clustering, to count the photoelectrons
                "time_resolution_in_ns": 0,
                "scintillators": sources,
            }
        },
    )
    out = _run_opt(
        monkeypatch,
        config_path,
        two_sources_stp,
        {"lar": optmap, "pen": optmap},
        tmp_path / "opt.lh5",
    )

    n_pes = ak.to_numpy(ak.num(out.dt))
    assert np.all(n_pes <= 40)
    assert out.is_saturated.to_list() == (n_pes == 40).tolist()
    assert 0 < np.sum(out.is_saturated) < len(out)


def _uniform_light(xloc, yloc, zloc, *, prob):
    # the steps of two_sources_stp are at (0.1, 0.1, 0.1) m
    assert np.allclose(np.stack([xloc, yloc, zloc]), 100)
    return np.full_like(xloc, prob)


def test_opt_script_light_model(
    tmp_path, monkeypatch, legend_testdata, l200_config_factory, two_sources_stp
):
    optmap = legend_testdata.get_path("remage/l200cfg01-optmap-dummy.lh5")
    config_path = l200_config_factory(
        tmp_path,
        {
            "opt": {
                "store_expected_pes": True,
                "max_pes_per_hit_combined": 10_000,
                "scintillators": {
                    "lar": {"table": "liquid_argon", "material": "lar"},
                    "pen": {
                        "table": "pen",
                        "material": "pen",
                        # does not change the expectation
                        "optmap_scaling_factor": 0.5,
                        "model": {
                            "function": f"{__name__}._uniform_light",
                            "kwargs": {"prob": 0.01},
                        },
                    },
                },
            }
        },
    )

    # the PEN source needs no map
    with pytest.raises(Exception, match="one optical map per light source"):
        _run_opt(
            monkeypatch,
            config_path,
            two_sources_stp,
            {"lar": optmap, "pen": optmap},
            tmp_path / "too_many_maps.lh5",
        )

    out = _run_opt(
        monkeypatch, config_path, two_sources_stp, {"lar": optmap}, tmp_path / "opt.lh5"
    )
    assert out.evtid.to_list() == [0, 1, 2, 3, 3]
    # PEN hits have three steps of 100 keV each, about 1600 photons, so 16
    # photoelectrons. Event 2 has two PEN hits, the second window of event 3 one
    expected = ak.to_numpy(out.expected_pes[[2, 4]])
    assert np.all(expected > 0)
    assert expected[0] / expected[1] == pytest.approx(2, rel=0.1)
    assert ak.all(ak.sum(out.energy[[2, 4]], axis=-1) > 0)
    assert ak.all(out.dt[:, 1:] >= out.dt[:, :-1])


def test_opt_script_light_model_per_sipm(
    tmp_path,
    monkeypatch,
    caplog,
    legend_testdata,
    l200_config_factory,
    two_sources_stp,
):
    optmap = legend_testdata.get_path("remage/l200cfg01-optmap-dummy.lh5")
    pen = {
        "table": "pen",
        "material": "pen",
        "model": {"function": f"{__name__}._uniform_light", "kwargs": {"prob": 0.01}},
    }
    config_path = l200_config_factory(
        tmp_path / "both",
        {
            "opt": {
                "scintillators": {
                    "lar": {"table": "liquid_argon", "material": "lar"},
                    "pen": pen,
                }
            }
        },
    )
    with caplog.at_level(logging.WARNING):
        out = _run_opt(
            monkeypatch,
            config_path,
            two_sources_stp,
            {"lar": optmap},
            tmp_path / "opt.lh5",
            per_sipm=True,
        )
    assert "['pen'] have a light model" in caplog.text
    # the windows with argon hits only
    assert out.evtid.to_list() == [0, 1, 3]

    config_path = l200_config_factory(
        tmp_path / "pen", {"opt": {"scintillators": {"pen": pen}}}
    )
    with pytest.raises(SimflowConfigError, match="all light sources"):
        _run_opt(
            monkeypatch,
            config_path,
            two_sources_stp,
            {},
            tmp_path / "pen.lh5",
            per_sipm=True,
        )
