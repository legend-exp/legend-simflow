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


@pytest.fixture(scope="session", autouse=True)
def dummyprod_geom_special_metadata():
    """Inject a writable ``special_metadata`` into each public-geometry config.

    All ``tests/dummyprod/inputs/simprod/config/geom/*-geom-config.yaml`` set
    ``public_geom: true``, which triggers
    ``pygeoml200.metadata.PublicMetadataProxy.update_special_metadata`` to
    mutate the dict returned by ``pygeoml200.core.configs.on(...)``. Since
    ``dbetto >=1.3.5`` returns ``TextDB.on(...)`` as read-only, this raises.
    Pre-resolve the default extra-metadata bundle to a sibling YAML and point
    each geom config at it via the ``special_metadata`` key, so
    ``pygeomtools.utils.load_dict_from_config`` wraps it in a writable
    ``AttrsDict``. Originals are restored on teardown.
    """
    geom_dir = testprod / "inputs/simprod/config/geom"
    extra_meta_path = geom_dir / "special_metadata.yaml"
    timestamp = "20230311T235840Z"
    extra_meta_path.write_text(
        yaml.safe_dump(core.configs.on(timestamp).to_dict(), sort_keys=False)
    )
    geom_configs = list(geom_dir.glob("*-geom-config.yaml"))
    backups = {p: p.read_text() for p in geom_configs}
    for p in geom_configs:
        cfg = yaml.safe_load(backups[p]) or {}
        cfg["special_metadata"] = str(extra_meta_path.resolve())
        p.write_text(yaml.safe_dump(cfg, sort_keys=False))
    yield
    for p, original in backups.items():
        p.write_text(original)
    extra_meta_path.unlink(missing_ok=True)


@pytest.fixture(scope="session")
def test_generate_gdml(config):
    geom_config = config.metadata.simprod.config.geom["legend-geom-config"].to_dict()
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
