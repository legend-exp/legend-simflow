from __future__ import annotations

import logging
import shutil
from pathlib import Path

import dbetto
import pytest
from dbetto import AttrsDict, TextDB

from legendsimflow import aggregate, utils, validate
from legendsimflow.exceptions import SimflowConfigError


@pytest.fixture
def l200data(tmp_path, test_l200data):
    prod = tmp_path / "prod"
    for d in [
        "inputs/datasets/statuses",
        "inputs/hardware/configuration/channelmaps",
        "inputs/hardware/detectors",
        "inputs/dataprod/config/tier/dsp",
        "generated/tier/raw",
        "generated/tier/hit",
        "generated/tier/pet",
        "generated/par/hit/cal",
    ]:
        (prod / d).mkdir(parents=True)

    dbetto.utils.write_dict({}, prod / "inputs/datasets/runinfo.yaml")
    dbetto.utils.write_dict(
        {
            "paths": {
                f"{kind}_{tier}": f"$_/generated/{kind}/{tier}"
                for kind, tier in [
                    ("tier", "raw"),
                    ("tier", "hit"),
                    ("tier", "pet"),
                    ("par", "hit"),
                ]
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
    return prod


def test_check_spec(tmp_path):
    spec = {
        "a|b": {"c": str, "d": float},
        "e/f": Path,
        "g": {},
    }
    obj = {"b": {"c": "x", "d": 1}, "e": {"f": tmp_path}, "g": {"h": 1}}
    assert validate.check_spec(obj, spec) == []

    obj = {"a": {"c": 1}, "e": {"f": tmp_path / "nope"}, "g": 1}
    assert validate.check_spec(obj, spec, "top") == [
        "no valid alternative among 'a|b' in 'top':",
        "  top/a/c: must be of type str, found int",
        "  top/a/d: not found",
        "  top/b: not found",
        f"top/e/f: {tmp_path / 'nope'} does not exist",
        "top/g: must be a dictionary or a folder",
    ]

    (tmp_path / "x/y").mkdir(parents=True)
    dbetto.utils.write_dict({"k": 1.0}, tmp_path / "x/z.yaml")
    db = TextDB(tmp_path, lazy=True)
    assert validate.check_spec(db, {"x/y": {}, "x": {"z": {"k": float}}}) == []
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

    config["allow_invalid_l200data"] = True
    with caplog.at_level(logging.WARNING):
        validate.l200data_layout(config)
    assert "channelmaps: not found" in caplog.text


def test_l200data_layout_not_a_folder(fresh_config, tmp_path):
    fresh_config.paths["l200data"] = tmp_path / "nope"
    with pytest.raises(SimflowConfigError, match="is not a folder"):
        validate.l200data_layout(fresh_config)


def test_l200data_layout_skipped_without_l200data(fresh_config):
    fresh_config.paths.pop("l200data", None)
    validate.l200data_layout(fresh_config)


def test_detinfo_against_l200data(fresh_config, monkeypatch):
    config = fresh_config
    sim = {"usability": {"l200-p02-r000-phy": {"V00048A": "on", "V00048B": "off"}}}
    data = {"l200-p02-r000-phy": {"V00048A": {"usability": "on"}}}
    data["l200-p02-r000-phy"]["V00048B"] = {"usability": "off"}

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

    config["allow_invalid_l200data"] = True
    validate.detinfo_against_l200data(config, sim)


def test_detinfo_skipped_without_l200data(config):
    config.paths.pop("l200data", None)
    validate.detinfo_against_l200data(config, {})
