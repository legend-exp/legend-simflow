from __future__ import annotations

import logging
import shutil
from pathlib import Path

import dbetto
import h5py
import pytest
from dbetto import TextDB

from legendsimflow import aggregate, utils, validate
from legendsimflow.exceptions import SimflowConfigError
from legendsimflow.metadata import parse_runid, reference_cal_run, runinfo


def _spec_paths(spec, prefix=""):
    for key, sub in spec.items():
        if sub:
            yield from _spec_paths(sub, f"{prefix}{key}/")
        else:
            yield f"{prefix}{key}"


@pytest.fixture
def l200data(tmp_path, test_l200data, fresh_config):
    prod = tmp_path / "prod"
    for d in ["datasets", "hardware"]:
        shutil.copytree(fresh_config.paths.metadata / d, prod / "inputs" / d)
    for d in ["inputs/dataprod/config/tier/dsp", "generated/par/hit/cal"]:
        (prod / d).mkdir(parents=True)

    dbetto.utils.write_dict(
        {
            "paths": {
                "tier_raw": "$_/generated/tier/raw",
                "tier_hit": "$_/generated/tier/hit",
                "tier_pet": "$_/generated/tier/pet",
                "par_hit": "$_/generated/par/hit",
            }
        },
        prod / "dataflow-config.yaml",
    )
    dbetto.utils.write_dict(
        [{"valid_from": "20220101T000000Z", "apply": ["cal/pars.yaml"]}],
        prod / "generated/par/hit/validity.yaml",
    )
    dbetto.utils.write_dict({}, prod / "generated/par/hit/cal/pars.yaml")
    shutil.copy(
        utils.lookup_dsp_config(test_l200data / "v3.0.0"),
        prod / "inputs/dataprod/config/tier/dsp",
    )

    # LH5 files with all the objects in LH5_SPECS, for all HPGe detectors
    meta = fresh_config.metadata
    runids = aggregate.gen_list_of_all_runids(fresh_config)
    objects = {
        path for spec in validate.LH5_SPECS.values() for path in _spec_paths(spec)
    }
    chmap = meta.channelmap(runinfo(meta, min(runids)).start_key)
    rawids = [c.daq.rawid for c in chmap.values() if c.get("system") == "geds"]
    for runid in runids | {reference_cal_run(meta, r) for r in runids}:
        _, period, run, datatype = parse_runid(runid)
        for tier in ["raw", "hit", "pet"]:
            folder = prod / f"generated/tier/{tier}/{datatype}/p{period:02d}/r{run:03d}"
            folder.mkdir(parents=True)
            with h5py.File(folder / f"{runid}-tier_{tier}.lh5", "w") as f:
                for rawid in rawids:
                    for obj in objects:
                        path = obj.format(rawid=rawid)
                        if path not in f:
                            f.create_dataset(path, data=[0])
    return prod


def test_check_spec(tmp_path):
    spec = {
        "a|b": {"c": {}, "d": {}},
        "e/f": Path,
        "g": {"h": {}},
    }
    obj = {"b": {"c": "x", "d": 1}, "e": {"f": tmp_path}, "g": {"h": 1}}
    assert validate.check_spec(obj, spec) == []

    obj = {"a": {"c": 1}, "e": {"f": tmp_path / "nope"}, "g": 1}
    assert validate.check_spec(obj, spec, "top") == [
        "no valid alternative among 'a|b' in 'top':",
        "  top/a/d: not found",
        "  top/b: not found",
        f"top/e/f: {tmp_path / 'nope'} does not exist",
        "top/g: must be a dictionary or a folder",
    ]

    (tmp_path / "x/y").mkdir(parents=True)
    dbetto.utils.write_dict({"k": 1.0}, tmp_path / "x/z.yaml")
    db = TextDB(tmp_path, lazy=True)
    assert validate.check_spec(db, {"x/y": {}, "x": {"z": {"k": {}}}}) == []
    assert validate.check_spec(db, {"x/w": {}}) == ["x/w: not found"]


def test_l200data_layout(fresh_config, l200data, caplog):
    config = fresh_config
    config.paths["l200data"] = l200data
    validate.l200data_layout(config)

    shutil.rmtree(l200data / "inputs/hardware/configuration/channelmaps")
    shutil.rmtree(l200data / "generated/tier/pet")
    dbetto.utils.write_dict(
        [{"valid_from": "20220301T000000Z", "apply": ["cal/pars.yaml"]}],
        l200data / "generated/par/hit/validity.yaml",
    )
    with pytest.raises(SimflowConfigError) as err:
        validate.l200data_layout(config)
    msg = str(err.value)
    assert "hardware/configuration/channelmaps: not found" in msg
    assert "dataflow config/paths/tier_pet:" in msg
    # the parameter files are checked only if the dataflow config is valid
    assert "l200-p02-r001-phy" not in msg

    shutil.copytree(l200data / "generated/tier/hit", l200data / "generated/tier/pet")
    with pytest.raises(SimflowConfigError) as err:
        validate.l200data_layout(config)
    msg = str(err.value)
    assert "l200-p02-r001-phy: no parameter file" in msg
    assert "l200-p02-r002-phy" not in msg

    config["validate_l200data"] = False
    with caplog.at_level(logging.WARNING):
        validate.l200data_layout(config)
    assert "is not validated" in caplog.text


def test_l200data_layout_lh5(fresh_config, l200data):
    config = fresh_config
    config.paths["l200data"] = l200data
    tier = l200data / "generated/tier"

    with h5py.File(next((tier / "pet/phy/p02/r000").glob("*.lh5")), "a") as f:
        del f["evt/spms/rawid"]
    # only the first file of the first run is opened
    with h5py.File(next((tier / "pet/phy/p02/r001").glob("*.lh5")), "a") as f:
        del f["evt/spms/t0"]
    for file in (tier / "hit/phy/p02/r003").glob("*.lh5"):
        file.unlink()

    with pytest.raises(SimflowConfigError) as err:
        validate.l200data_layout(config)
    msg = str(err.value)
    assert "evt/spms/rawid: not found" in msg
    assert "spms/t0" not in msg
    assert f"l200-p02-r003-phy: no LH5 files in {tier / 'hit/phy/p02/r003'}" in msg


def test_l200data_layout_superpulses_t0(fresh_config, l200data, monkeypatch):
    config = fresh_config
    config.paths["l200data"] = l200data
    get_tier_settings = validate.get_tier_settings
    monkeypatch.setattr(
        validate,
        "get_tier_settings",
        lambda c, t: get_tier_settings(c, t)
        | ({"simulate_psd_with_psl": True} if t == "hit" else {}),
    )
    monkeypatch.setattr(validate.patterns, "compute_superpulses", lambda *_, **__: True)
    validate.l200data_layout(config)

    meta = config.metadata
    runs = {
        reference_cal_run(meta, r) for r in aggregate.gen_list_of_all_runids(config)
    }
    _, period, run, datatype = parse_runid(max(runs))
    folder = l200data / f"generated/tier/pet/{datatype}/p{period:02d}/r{run:03d}"
    file = next(folder.glob("*.lh5"))
    with h5py.File(file, "a") as f:
        del f["evt/spms/event_t0"]

    with pytest.raises(SimflowConfigError, match="no evt/spms/event_t0 or evt/geds/t0"):
        validate.l200data_layout(config)

    with h5py.File(file, "a") as f:
        f.create_dataset("evt/geds/t0", data=[0])
    validate.l200data_layout(config)


def test_l200data_layout_not_a_folder(fresh_config, tmp_path):
    fresh_config.paths["l200data"] = tmp_path / "nope"
    with pytest.raises(SimflowConfigError, match="is not a folder"):
        validate.l200data_layout(fresh_config)


def test_l200data_layout_skipped_without_l200data(fresh_config):
    fresh_config.paths.pop("l200data", None)
    validate.l200data_layout(fresh_config)


def test_l200data_layout_usabilities(fresh_config, l200data):
    config = fresh_config
    config.paths["l200data"] = l200data
    statuses = l200data / "inputs/datasets/statuses/l200-p02-r006-T%-all-config.yaml"
    status = dbetto.utils.load_dict(statuses)
    status["V99000A"]["usability"] = "off"
    status["V99000A"]["psd"]["status"]["low_aoe"] = "present"
    dbetto.utils.write_dict(status, statuses)

    with pytest.raises(SimflowConfigError) as err:
        validate.l200data_layout(config)
    msg = str(err.value)
    assert (
        "l200-p02-r001-phy V99000A: usability is 'on' in the Simflow metadata, "
        "'off' in the data production metadata"
    ) in msg
    assert (
        "l200-p02-r001-phy V99000A: psd_usability is 'valid' in the Simflow "
        "metadata, 'present' in the data production metadata"
    ) in msg


def test_input_paths(fresh_config, tmp_path):
    config = fresh_config
    validate.input_paths(config)

    config.paths["config"] = tmp_path / "nope"
    config.paths["hpge_impurity_model"] = tmp_path / "impurity.yaml"
    with pytest.raises(SimflowConfigError) as e:
        validate.input_paths(config)
    assert "paths/config" in str(e.value)
    assert "paths/hpge_impurity_model" in str(e.value)


def test_input_paths_optical_maps(fresh_config, tmp_path):
    config = fresh_config
    config.paths.optical_maps["lar"] = tmp_path / "lar.lh5"
    config.paths.optical_maps["lar_patch"] = {"sim1": tmp_path / "patch.lh5"}

    # the maps are not read if the opt tier is not built
    config["make_steps"] = ["vtx", "stp"]
    validate.input_paths(config)

    config["make_steps"] = ["vtx", "stp", "opt"]
    with pytest.raises(SimflowConfigError) as e:
        validate.input_paths(config)
    assert "paths/optical_maps/lar:" in str(e.value)
    assert "paths/optical_maps/lar_patch/sim1" in str(e.value)

    (tmp_path / "lar.lh5").touch()
    (tmp_path / "patch.lh5").touch()
    validate.input_paths(config)


@pytest.mark.parametrize("system", ["geds", "spms"])
def test_l200data_layout_no_detector_on(fresh_config, l200data, monkeypatch, system):
    config = fresh_config
    config.paths["l200data"] = l200data
    monkeypatch.setattr(validate, "_usability_problems", lambda *_: [])
    chmap = config.metadata.channelmap
    monkeypatch.setattr(
        config.metadata,
        "channelmap",
        lambda *a, **kw: {
            k: v for k, v in chmap(*a, **kw).items() if v.get("system") != system
        },
    )
    if system == "geds":
        # the HPGe raw files would be read if an HPGe detector was on
        for file in (l200data / "generated/tier/raw").rglob("*.lh5"):
            file.unlink()
    validate.l200data_layout(config)
