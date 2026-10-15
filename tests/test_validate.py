from __future__ import annotations

from pathlib import Path

import pytest
from dbetto import AttrsDict

from legendsimflow import aggregate, validate
from legendsimflow.exceptions import SimflowConfigError


def test_detinfo_against_l200data(config, monkeypatch):
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

    config["allow_metadata_differences"] = True
    validate.detinfo_against_l200data(config, sim)


def test_detinfo_skipped_without_l200data(config):
    config.paths.pop("l200data", None)
    validate.detinfo_against_l200data(config, {})
