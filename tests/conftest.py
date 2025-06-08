from pathlib import Path

import legenddataflowscripts
import pytest
import yaml

testprod = Path(__file__).parent / "test-prod"
config_filename = testprod / "simflow-config.yaml"


@pytest.fixture(scope="session")
def simflow_config():
    with config_filename.open() as f:
        return legenddataflowscripts.workflow.utils.subst_vars(
            yaml.safe_load(f), var_values={"_": str(Path(__file__).parent)}
        )


class mock_workflow_class:
    def __init__(self):
        self.overwrite_configfiles = [config_filename]


@pytest.fixture(scope="module")
def mock_workflow():
    return mock_workflow_class()


# @pytest.fixture
# def mock_os_environ(monkeypatch):
#     monkeypatch.setenv("PRODENV", "prod")
#     return os.environ
