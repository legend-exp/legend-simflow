from __future__ import annotations

import contextlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

import dbetto.utils
import pytest
import yaml

from legendsimflow import aggregate, utils
from legendsimflow.scripts import extract_hpge_current_pulse_model
from legendsimflow.scripts.make_simstat_partition_file import main as simstat_main
from legendsimflow.scripts.pars import extract_hpge_observables_models

testprod = Path(__file__).parent.parent / "dummyprod"
repo_root = Path(__file__).parent.parent.parent

_RUNIDS_L1000 = ("l200-p03-r000-phy", "l200-p03-r001-phy")


@contextlib.contextmanager
def _override_argv(*args):
    """Temporarily replace sys.argv for session-scoped fixture invocations."""
    old = sys.argv[:]
    sys.argv = list(args)
    try:
        yield
    finally:
        sys.argv = old


@pytest.fixture(scope="session")
def legend_gdml_path(tmp_path_factory):
    """Generate the legend GDML file using legend-pygeom-l200.

    Calls the ``legend-pygeom-l200`` CLI with the dummyprod geometry config and
    metadata, writing a pygeomtools-compatible GDML to a session-scoped
    temporary directory. The result is cached for the full test session (~5 s
    one-time cost).
    """
    out_dir = tmp_path_factory.mktemp("legend_gdml")
    gdml_path = out_dir / "legend.gdml"
    geom_config = testprod / "inputs/simprod/config/geom/legend-geom-config.yaml"
    env = os.environ.copy()
    env["LEGEND_METADATA"] = str(testprod / "inputs")
    subprocess.run(
        [
            "legend-pygeom-l200",
            "--config",
            str(geom_config),
            "--",
            str(gdml_path),
        ],
        check=True,
        env=env,
    )
    return gdml_path


@pytest.fixture(scope="session")
def legend_stp_path(tmp_path_factory, legend_gdml_path):
    """Generate a legend stp LH5 file by running remage with the public geometry.

    Skips if remage is not installed (requires the pixi test environment).  The
    simulation is a minimal 2 MeV gamma source confined to all V-type detectors,
    producing hits in several germanium detectors.  The result is cached for the
    full test session.
    """
    if shutil.which("remage") is None:
        pytest.skip("remage not installed")

    out_dir = tmp_path_factory.mktemp("legend_stp")
    stp_file = out_dir / "legend-test_hit_sim-job_0000-tier_stp.lh5"

    # Each macro command is a separate list element, matching how commands.py
    # builds the remage invocation with macro_free=True.
    commands = [
        "/RMG/Manager/Randomization/Seed 42",
        "/RMG/Geometry/RegisterDetectorsFromGDML Germanium",
        "/RMG/Geometry/RegisterDetectorsFromGDML Scintillator",
        "/RMG/Geometry/GDMLDisableOverlapCheck",
        "/RMG/Output/NtupleUseVolumeName true",
        "/run/initialize",
        "/RMG/Output/Germanium/StoreSinglePrecisionPosition",
        "/RMG/Output/Germanium/StoreSinglePrecisionEnergy",
        "/RMG/Output/Vertex/StoreSinglePrecisionPosition",
        "/RMG/Output/Germanium/StoreTrackID true",
        "/RMG/Output/Germanium/EdepCutLow 1 eV",
        "/RMG/Generator/Select GPS",
        "/gps/particle gamma",
        "/gps/energy 2 MeV",
        "/gps/ang/type iso",
        "/RMG/Generator/Confine Volume",
        "/RMG/Generator/Confinement/Physical/AddVolume V.*",
        "/RMG/Generator/Confinement/MaxSamplingTrials 100000",
        "/run/beamOn 10000",
    ]

    subprocess.run(
        [
            "remage",
            "--gdml-files",
            str(legend_gdml_path),
            "--output-file",
            str(stp_file),
            "--",
            *commands,
        ],
        check=True,
    )

    return stp_file


@pytest.fixture(scope="session")
def legend_dtmap_path(tmp_path_factory):
    """Generate a drift time map LH5 file for V05261B using l1000dsg01 metadata.

    Calls the ``make_hpge_drift_time_maps.jl`` Julia script with the l1000dsg01
    dtmap settings (coarse 10 mm grid for speed).  Skips if Julia is not
    installed.  The result is cached for the full test session and can be used
    as ``--dtmap-files`` input to the hit tier script.
    """
    if shutil.which("julia") is None:
        pytest.skip("julia not installed")

    out_dir = tmp_path_factory.mktemp("legend_dtmap")
    dtmap_file = out_dir / "V05261B-4200V-hpge-drift-time-map.lh5"

    julia_script = (
        repo_root / "workflow/src/legendsimflow/scripts/make_hpge_drift_time_maps.jl"
    )
    julia_project = repo_root / "workflow/src/LegendSimflow.jl"
    dtmap_settings = (
        testprod / "inputs/simprod/config/pars/l1000dsg01/geds/dtmap/settings.yaml"
    )

    subprocess.run(
        [
            "julia",
            "--project=" + str(julia_project),
            "--threads",
            "1",
            str(julia_script),
            "--detector",
            "V05261B",
            "--metadata",
            str(testprod / "inputs"),
            "--opv",
            "4200",
            "--dtmap-settings",
            str(dtmap_settings),
            "--output-file",
            str(dtmap_file),
        ],
        check=True,
        cwd=repo_root,
    )
    return dtmap_file


def _l1000_config(tmp_dir: Path) -> Path:
    """Write a minimal simflow-config-l1000.yaml to *tmp_dir* and return its path.

    Only ``paths.metadata`` is overridden to point at the dummyprod inputs;
    all other path entries use ``$_`` substitution resolved to *tmp_dir*.
    """
    raw = yaml.safe_load((testprod / "simflow-config-l1000.yaml").read_text())
    raw["paths"]["metadata"] = str(testprod / "inputs")
    config_path = tmp_dir / "simflow-config-l1000.yaml"
    config_path.write_text(yaml.safe_dump(raw))
    return config_path


@pytest.fixture(scope="session")
def legend_currmod_paths(tmp_path_factory):
    """Run ``extract_hpge_current_pulse_model`` for both l1000dsg01 runids.

    The l1000dsg01 metadata has a ``default`` key so no l200data is required.
    Returns a dict mapping each runid to its output YAML path.
    """
    out_dir = tmp_path_factory.mktemp("legend_currmod")
    config_path = _l1000_config(out_dir)

    paths = {}
    for runid in _RUNIDS_L1000:
        pars_file = out_dir / f"{runid}-model.yaml"
        plot_file = out_dir / f"{runid}-fit-results.pdf"

        with _override_argv(
            "extract-hpge-currmod",
            "--runid",
            runid,
            "--pars-file",
            str(pars_file),
            "--plot-file",
            str(plot_file),
            "--simflow-config",
            str(config_path),
        ):
            extract_hpge_current_pulse_model.main()

        paths[runid] = pars_file

    return paths


@pytest.fixture(scope="session")
def legend_simstat_part_path(tmp_path_factory, legend_stp_path):
    """Run ``make_simstat_partition_file`` with both l1000dsg01 runids.

    Depends on ``legend_stp_path``; skips if remage is not installed.
    Returns the path to the output partition YAML file.
    """
    out_dir = tmp_path_factory.mktemp("legend_simstat")
    config_path = _l1000_config(out_dir)
    output_file = out_dir / "partitions.yaml"

    with _override_argv(
        "make-simstat-partition-file",
        "--stp-files",
        str(legend_stp_path),
        "--runlist",
        *_RUNIDS_L1000,
        "--output-file",
        str(output_file),
        "--simflow-config",
        str(config_path),
    ):
        simstat_main()

    return output_file


@pytest.fixture(scope="session")
def legend_detector_usabilities_path(tmp_path_factory):
    """Cache detector usabilities for all l1000dsg01 runs.

    Calls ``aggregate.gen_list_of_all_usabilities`` with the l1000dsg01
    metadata (no remage or Julia required).  Returns the path to the output
    YAML file.
    """
    out_dir = tmp_path_factory.mktemp("legend_usabilities")
    config_path = _l1000_config(out_dir)
    output_file = out_dir / "detector_usabilities.yaml"

    config = utils.init_simflow_context(config_path, workflow=None).config
    dbetto.utils.write_dict(
        aggregate.gen_list_of_all_usabilities(config).to_dict(), output_file
    )

    return output_file


@pytest.fixture(scope="session")
def legend_hpge_obs_paths(tmp_path_factory):
    """Run ``extract_hpge_observables_models`` for both l1000dsg01 runids.

    The l1000dsg01 metadata has ``default`` keys for all three observables so
    no l200data is required.  Returns a dict
    ``{runid: {eresmod: path, aoeresmod: path, psdcuts: path}}``.
    """
    out_dir = tmp_path_factory.mktemp("legend_hpge_obs")
    config_path = _l1000_config(out_dir)

    paths = {}
    for runid in _RUNIDS_L1000:
        eresmod_file = out_dir / f"{runid}-eresmod.yaml"
        aoeresmod_file = out_dir / f"{runid}-aoeresmod.yaml"
        psdcuts_file = out_dir / f"{runid}-psdcuts.yaml"

        with _override_argv(
            "extract-hpge-obs-models",
            "--runid",
            runid,
            "--eresmod-file",
            str(eresmod_file),
            "--aoeresmod-file",
            str(aoeresmod_file),
            "--psdcuts-file",
            str(psdcuts_file),
            "--simflow-config",
            str(config_path),
        ):
            extract_hpge_observables_models.main()

        paths[runid] = {
            "eresmod": eresmod_file,
            "aoeresmod": aoeresmod_file,
            "psdcuts": psdcuts_file,
        }

    return paths
