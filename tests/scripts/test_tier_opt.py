from __future__ import annotations

import logging
import sys
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

from legendsimflow.exceptions import SimflowConfigError
from legendsimflow.scripts.tier import opt

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
        "--optmap-lar", str(optmap_path),
        "--geom-file", str(legend_gdml_path),
        "--usability-file", str(legend_detector_usabilities_path / "usability.yaml"),
        "--jobid", "0000",
        "--scintillator-volume-name", "liquid_argon",
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
        "--scintillator-volume-name", "liquid_argon",
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
        "optmap-lar": "optmap.lh5",
        "geom-file": "geom.gdml",
        "simstat-part-file": "partitions.yaml",
    } | {k.replace("_", "-"): v for k, v in paths.items()}
    return [
        "opt",
        *[a for k, v in files.items() for a in (f"--{k}", str(v))],
        "--usability-file", str(usability_dir / "usability.yaml"),
        "--jobid", "0000",
        "--scintillator-volume-name", "liquid_argon",
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
    i = argv.index("--optmap-lar")
    del argv[i : i + 2]
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SimflowConfigError, match="--optmap-lar is required"):
        opt.main()


@pytest.mark.needs_remage
def test_opt_script_unknown_scintillator_volume(
    tmp_path, monkeypatch, legend_gdml_path, legend_detector_usabilities_path
):
    argv = _argv(
        tmp_path, legend_detector_usabilities_path, geom_file=str(legend_gdml_path)
    )
    argv[argv.index("--scintillator-volume-name") + 1] = "not_a_volume"
    monkeypatch.setattr(sys, "argv", argv)
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
        optmap_lar=legend_testdata.get_path("remage/l200cfg01-optmap-dummy.lh5"),
        geom_file=str(legend_gdml_path),
        simstat_part_file=str(part_file),
    ))  # fmt: skip

    # nothing is written, so the TCM cannot be built at the end
    with caplog.at_level(logging.WARNING), pytest.raises(FileNotFoundError):
        opt.main()
    assert not (tmp_path / "opt.lh5").exists()
    assert any(
        "liquid_argon not found" in r.message and r.levelno == logging.WARNING
        for r in caplog.records
    )
