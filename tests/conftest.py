from __future__ import annotations

import shutil
from pathlib import Path

import legenddataflowscripts
import pytest
import yaml
from dbetto import AttrsDict
from legendmeta import LegendMetadata
from legendtestdata import LegendTestData
from pygeoml200 import core

from legendsimflow.utils import apply_path_defaults

testprod = Path(__file__).parent / "dummyprod"
config_filename = testprod / "simflow-config.yaml"


@pytest.fixture(scope="session")
def legend_testdata():
    ldata = LegendTestData()
    ldata.checkout("68d8b49")
    return ldata


@pytest.fixture(scope="session", autouse=True)
def dummyprod_optmap(legend_testdata):
    """Copy the real optical map from legend_testdata into dummyprod.

    All configs referencing ``$_/inputs/simprod/l200cfg01-optmap-dummy.lh5``
    find a valid LH5 file.  The file is gitignored; this fixture is the sole
    source of truth for test runs.
    """
    src = Path(legend_testdata.get_path("remage/l200cfg01-optmap-dummy.lh5"))
    dst = testprod / "inputs/simprod/l200cfg01-optmap-dummy.lh5"
    shutil.copy2(src, dst)
    yield
    dst.unlink(missing_ok=True)


@pytest.fixture(scope="session")
def test_generate_gdml(config):
    geom_config = config.metadata.simprod.config.geom["legend-geom-config"]

    return core.construct(
        use_detailed_fiber_model=False, config=geom_config, public_geometry=True
    )


def make_config():
    with config_filename.open() as f:
        config = yaml.safe_load(f)

    legenddataflowscripts.subst_vars(config, var_values={"_": testprod})
    assert config is not None

    # convert all strings in the "paths" block to pathlib.Path
    def _make_path(d):
        for k, v in d.items():
            if isinstance(v, str):
                d[k] = Path(v)
            else:
                d[k] = _make_path(v)
        return d

    config["paths"] = _make_path(config["paths"])
    apply_path_defaults(config["paths"])

    metadata = LegendMetadata(testprod / "inputs")

    config["metadata"] = metadata
    config["_proctime"] = "now"

    return AttrsDict(config)


@pytest.fixture(scope="session")
def config():
    return make_config()


@pytest.fixture
def fresh_config():
    return make_config()


class mock_workflow_class:
    def __init__(self):
        self.overwrite_configfiles = [config_filename]


@pytest.fixture(scope="module")
def mock_workflow():
    return mock_workflow_class()


@pytest.fixture(scope="session")
def test_l200data():
    return Path(__file__).parent / "l200data"
