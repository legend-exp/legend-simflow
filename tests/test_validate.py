from __future__ import annotations

import logging
import shutil
from pathlib import Path

import dbetto
import h5py
import pytest
from dbetto import AttrsDict, TextDB

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
    for d in [
        "inputs/datasets/statuses",
        "inputs/hardware/configuration/channelmaps",
        "inputs/hardware/detectors",
        "inputs/dataprod/config/tier/dsp",
        "generated/par/hit/cal",
    ]:
        (prod / d).mkdir(parents=True)

    dbetto.utils.write_dict({}, prod / "inputs/datasets/runinfo.yaml")
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


def test_l200data_layout_not_a_folder(fresh_config, tmp_path):
    fresh_config.paths["l200data"] = tmp_path / "nope"
    with pytest.raises(SimflowConfigError, match="is not a folder"):
        validate.l200data_layout(fresh_config)


def test_l200data_layout_skipped_without_l200data(fresh_config):
    fresh_config.paths.pop("l200data", None)
    validate.l200data_layout(fresh_config)


def test_detinfo_against_l200data(fresh_config, monkeypatch, caplog):
    config = fresh_config
    sim = {"usability": {"l200-p02-r000-phy": {"V00048A": "on", "V00048B": "off"}}}
    data = {
        "l200-p02-r000-phy": {
            "V00048A": {"usability": "on"},
            "V00048B": {"usability": "off"},
        }
    }

    monkeypatch.setattr(validate, "LegendMetadata", lambda *_, **__: None)
    queried = []

    def fake_usabilities(_, runids):
        queried.append(set(runids))
        return AttrsDict(data)

    monkeypatch.setattr(aggregate, "gen_list_of_all_usabilities", fake_usabilities)

    config.paths["l200data"] = Path("/does/not/matter")
    validate.detinfo_against_l200data(config, sim)
    assert queried == [{"l200-p02-r000-phy"}]

    data["l200-p02-r000-phy"]["V00048B"]["usability"] = "on"
    with pytest.raises(SimflowConfigError) as err:
        validate.detinfo_against_l200data(config, sim)
    assert (
        "l200-p02-r000-phy V00048B: Simflow metadata has 'off', "
        "data production has 'on'"
    ) in str(err.value)
    assert "V00048A" not in str(err.value)

    config["validate_l200data"] = False
    with caplog.at_level(logging.WARNING):
        validate.detinfo_against_l200data(config, sim)
    assert "is not validated" in caplog.text


def test_detinfo_skipped_without_l200data(config):
    config.paths.pop("l200data", None)
    validate.detinfo_against_l200data(config, {})
