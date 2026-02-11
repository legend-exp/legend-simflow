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
    ldata.checkout("0bb2dfc")
    return ldata


@pytest.fixture(scope="session")
def legend_test_metadata(legend_testdata):
    return LegendMetadata(legend_testdata["legend/metadata"])


@pytest.fixture(scope="session")
def test_generate_gdml(config):
    geom_config = config.metadata.simprod.config.geom["l200p02-geom-config"]

    return pygeoml200.core.construct(
        use_detailed_fiber_model=False, config=geom_config, public_geometry=True
    )


def make_config(legend_testdata):
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

    def _copy_skip_existing(src, dst):
        if not Path(dst).exists():
            shutil.copy2(src, dst)

    for fd in ("hardware", "datasets"):
        shutil.copytree(
            legend_testdata[f"legend/metadata/{fd}"],
            testprod / "inputs" / fd,
            copy_function=_copy_skip_existing,
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


@pytest.fixture(scope="session")
def test_l200data():
    return Path(__file__).parent / "dummyprod/l200data"
