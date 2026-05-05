from __future__ import annotations

import contextlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

import dbetto.utils
import h5py
import pytest
import yaml

from legendsimflow import aggregate, utils
from legendsimflow.scripts import extract_hpge_current_pulse_model
from legendsimflow.scripts.make_simstat_partition_file import main as simstat_main
from legendsimflow.scripts.pars import extract_hpge_observables_models
from legendsimflow.scripts.tier import cvt, evt, hit, opt

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
def legend_dtmap_path():
    """Return the path to the pre-built dummy drift time map for V05261B.

    The file lives in ``dummyprod/inputs/simprod/`` and contains constant
    1000 ns drift times on a 1 mm grid covering the V05261B detector volume.
    It is used as ``--dtmap-files`` input to the hit tier script.
    """
    return testprod / "inputs/simprod/V05261B-4200V-hpge-drift-time-map.lh5"


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
    For each runid, runs the per-detector extraction script over all modelable
    HPGes and merges the outputs into a single ``{runid}-model.yaml`` keyed by
    detector name (mirroring the ``merge_current_pulse_model_pars`` rule).
    Returns a dict mapping each runid to its merged YAML path.
    """
    out_dir = tmp_path_factory.mktemp("legend_currmod")
    config_path = _l1000_config(out_dir)
    config = utils.init_simflow_context(str(config_path), workflow=None).config

    paths = {}
    for runid in _RUNIDS_L1000:
        merged_file = out_dir / f"{runid}-model.yaml"
        merged: dict = {}
        for hpge in aggregate.gen_list_of_hpges_valid_for_modeling(config, runid):
            pars_file = out_dir / f"{runid}-{hpge}-model.yaml"
            plot_file = out_dir / f"{runid}-{hpge}-fit-result.pdf"
            with _override_argv(
                "extract-hpge-currmod",
                "--runid",
                runid,
                "--hpge-detector",
                hpge,
                "--pars-file",
                str(pars_file),
                "--plot-file",
                str(plot_file),
                "--simflow-config",
                str(config_path),
            ):
                extract_hpge_current_pulse_model.main()
            merged[hpge] = dbetto.utils.load_dict(pars_file)

        dbetto.utils.write_dict(merged, merged_file)
        paths[runid] = merged_file

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


@pytest.fixture(scope="session")
def legend_opt_path(
    tmp_path_factory,
    legend_testdata,
    legend_stp_path,
    legend_gdml_path,
    legend_simstat_part_path,
    legend_detector_usabilities_path,
):
    """Run ``opt.main()`` and return the path to the output opt LH5 file.

    Uses the dummy optical map from legend-testdata.  Skips if remage is not
    installed (the stp fixture already enforces this).
    """
    out_dir = tmp_path_factory.mktemp("legend_opt")
    config_path = _l1000_config(out_dir)
    opt_file = out_dir / "opt.lh5"

    optmap_path = Path(legend_testdata.get_path("remage/l200cfg01-optmap-dummy.lh5"))

    with _override_argv(
        "opt",
        "--stp-file",
        str(legend_stp_path),
        "--optmap-lar",
        str(optmap_path),
        "--geom-file",
        str(legend_gdml_path),
        "--simstat-part-file",
        str(legend_simstat_part_path),
        "--detector-usabilities-file",
        str(legend_detector_usabilities_path),
        "--jobid",
        "0000",
        "--opt-file",
        str(opt_file),
        "--scintillator-volume-name",
        "liquid_argon",
        "--simflow-config",
        str(config_path),
    ):
        opt.main()

    return opt_file


@pytest.fixture(scope="session")
def legend_hit_path(
    tmp_path_factory,
    legend_stp_path,
    legend_gdml_path,
    legend_dtmap_path,
    legend_currmod_paths,
    legend_hpge_obs_paths,
    legend_simstat_part_path,
    legend_detector_usabilities_path,
):
    """Run ``hit.main()`` and return the path to the output hit LH5 file.

    Sets up a full pars directory with energy-resolution models, A/E models,
    PSD cuts, current-pulse models, and drift-time maps for the two l1000dsg01
    run IDs.  Skips if remage or Julia are not installed (the stp / dtmap
    fixtures already enforce this).
    """
    out_dir = tmp_path_factory.mktemp("legend_hit")
    pars_dir = out_dir / "pars"

    for runid in _RUNIDS_L1000:
        obs = legend_hpge_obs_paths[runid]
        for subdir, src, dest_name in (
            ("hpge/eresmod", obs["eresmod"], f"{runid}-model.yaml"),
            ("hpge/aoeresmod", obs["aoeresmod"], f"{runid}-model.yaml"),
            ("hpge/psdcuts", obs["psdcuts"], f"{runid}-psd-cuts.yaml"),
            ("hpge/currmod", legend_currmod_paths[runid], f"{runid}-model.yaml"),
        ):
            dest_dir = pars_dir / subdir
            dest_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy(src, dest_dir / dest_name)

    dtmap_dir = pars_dir / "hpge/dtmaps"
    dtmap_dir.mkdir(parents=True)
    # r000: real dtmap for V05261B → PSD will be computed
    shutil.copy(
        legend_dtmap_path,
        dtmap_dir / f"{_RUNIDS_L1000[0]}-hpge-drift-time-maps.lh5",
    )
    # r001: empty dtmap → dt_map = None → PSD will be NaN
    with h5py.File(dtmap_dir / f"{_RUNIDS_L1000[1]}-hpge-drift-time-maps.lh5", "w"):
        pass

    raw = yaml.safe_load((testprod / "simflow-config-l1000.yaml").read_text())
    raw["paths"]["metadata"] = str(testprod / "inputs")
    raw["paths"]["pars"] = str(pars_dir)
    config_path = out_dir / "simflow-config-l1000.yaml"
    config_path.write_text(yaml.safe_dump(raw))

    hit_file = out_dir / "hit.lh5"

    with _override_argv(
        "hit",
        "--stp-file",
        str(legend_stp_path),
        "--jobid",
        "0000",
        "--hit-file",
        str(hit_file),
        "--geom-file",
        str(legend_gdml_path),
        "--dtmap-files",
        str(dtmap_dir / f"{_RUNIDS_L1000[0]}-hpge-drift-time-maps.lh5"),
        "--currmod-files",
        str(legend_currmod_paths[_RUNIDS_L1000[0]]),
        str(legend_currmod_paths[_RUNIDS_L1000[1]]),
        "--simstat-part-file",
        str(legend_simstat_part_path),
        "--detector-usabilities-file",
        str(legend_detector_usabilities_path),
        "--simflow-config",
        str(config_path),
    ):
        hit.main()

    return hit_file


@pytest.fixture(scope="session")
def legend_evt_path(
    tmp_path_factory,
    legend_stp_path,
    legend_opt_path,
    legend_hit_path,
    legend_simstat_part_path,
    legend_detector_usabilities_path,
):
    """Run ``evt.main()`` and return the path to the output evt LH5 file.

    Skips if remage is not installed (the stp fixture already enforces this).
    """
    out_dir = tmp_path_factory.mktemp("legend_evt")
    config_path = _l1000_config(out_dir)
    evt_file = out_dir / "evt.lh5"

    with _override_argv(
        "evt",
        "--stp-file",
        str(legend_stp_path),
        "--opt-file",
        str(legend_opt_path),
        "--hit-file",
        str(legend_hit_path),
        "--simstat-part-file",
        str(legend_simstat_part_path),
        "--detector-usabilities-file",
        str(legend_detector_usabilities_path),
        "--jobid",
        "0000",
        "--evt-file",
        str(evt_file),
        "--simflow-config",
        str(config_path),
    ):
        evt.main()

    return evt_file


@pytest.fixture(scope="session")
def legend_cvt_path(tmp_path_factory, legend_evt_path):
    """Run ``cvt.main()`` and return the path to the output cvt LH5 file.

    Skips if remage is not installed (the evt fixture already enforces this).
    """
    out_dir = tmp_path_factory.mktemp("legend_cvt")
    config_path = _l1000_config(out_dir)
    cvt_file = out_dir / "cvt.lh5"

    with _override_argv(
        "cvt",
        "--evt-files",
        str(legend_evt_path),
        "--cvt-file",
        str(cvt_file),
        "--simflow-config",
        str(config_path),
    ):
        cvt.main()

    return cvt_file
