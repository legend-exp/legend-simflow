from __future__ import annotations

import shutil
from pathlib import Path

import legenddataflowscripts
import pytest
import yaml
from dbetto import AttrsDict
from legendmeta import LegendMetadata
from legendtestdata import LegendTestData

testprod = Path(__file__).parent / "dummyprod"
config_filename = testprod / "simflow-config.yaml"


@pytest.fixture(scope="session")
def legend_testdata():
    ldata = LegendTestData()
    ldata.checkout("8247690")
    return ldata


@pytest.fixture(scope="session")
def legend_test_metadata(legend_testdata):
    return LegendMetadata(legend_testdata["legend/metadata"])


def make_config(legend_testdata):
    with config_filename.open() as f:
        config = yaml.safe_load(f)

    legenddataflowscripts.subst_vars(config, var_values={"_": testprod})
    assert config is not None

    for fd in ("hardware", "datasets"):
        shutil.copytree(
            legend_testdata[f"legend/metadata/{fd}"],
            testprod / "inputs" / fd,
            dirs_exist_ok=True,
        )

    metadata = LegendMetadata(testprod / "inputs")

    config["metadata"] = metadata

    return AttrsDict(config)


@pytest.fixture(scope="session")
def config(legend_testdata):
    return make_config(legend_testdata)


@pytest.fixture
def fresh_config(legend_testdata):
    return make_config(legend_testdata)


class mock_workflow_class:
    def __init__(self):
        self.overwrite_configfiles = [config_filename]


@pytest.fixture(scope="module")
def mock_workflow():
    return mock_workflow_class()
