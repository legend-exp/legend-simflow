from __future__ import annotations

import contextlib
import os
import shutil
import subprocess
import sys
from collections.abc import Mapping
from pathlib import Path

import dbetto
import dbetto.utils
import h5py
import lh5
import pytest
import yaml

from legendsimflow import aggregate, utils
from legendsimflow.scripts import extract_hpge_current_pulse_model
from legendsimflow.scripts.make_simstat_partition_file import main as simstat_main
from legendsimflow.scripts.pars import extract_hpge_observables_models
from legendsimflow.scripts.tier import cvt, evt, hit, opt

testprod = Path(__file__).parent.parent / "dummyprod"
repo_root = Path(__file__).parent.parent.parent

_RUNIDS_L200 = ("l200-p03-r000-phy", "l200-p03-r001-phy")


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
    """Build the array of the `l200cfg01` experiment, the one the full chain runs on."""
    out_dir = tmp_path_factory.mktemp("legend_gdml")
    gdml_path = out_dir / "legend.gdml"

    # read the configuration through the metadata database, which expands its
    # $_ paths, as the gen_geom_config rule does
    metadata = dbetto.TextDB(testprod / "legend-metadata", lazy=True)
    geom_config = out_dir / "geom-config.yaml"
    dbetto.utils.write_dict(
        metadata.simprod.config.geom["l200cfg01-geom-config"].to_dict(), geom_config
    )

    env = os.environ.copy()
    env["LEGEND_METADATA"] = str(testprod / "legend-metadata")
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
    """Simulate 2 MeV gammas in the germanium detectors of the array, with remage."""
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
def legend_dtmap_path(tmp_path_factory, legend_testdata):
    """Write the drift time maps of the two modelled detectors.

    legend-testdata ships the maps of `V99999Z`, of which every detector of the
    array is a copy, so the same maps serve all of them under their own name.
    """
    source = Path(
        legend_testdata.get_path("remage/V99999Z-3500V-hpge-drift-time-map.lh5")
    )
    maps = lh5.read("V99999Z", source)

    path = tmp_path_factory.mktemp("legend_dtmap") / "hpge-drift-time-maps.lh5"
    for detector in ("V00001A", "V00001B"):
        lh5.write(maps, detector, path, wo_mode="write_safe")

    return path


def _l200_config(tmp_dir: Path, settings_by_tier: Mapping | None = None) -> Path:
    """Write a minimal simflow-config-l200.yaml to *tmp_dir* and return its path.

    The metadata is the one the full-chain workflow test runs on, assembled by
    the ``dummyprod_testdata`` fixture.

    ``settings_by_tier`` maps a tier name to the settings keys to overwrite. When
    it is given, the metadata is copied to `<tmp_dir>/legend-metadata` first, so
    the committed tree stays untouched. That copy is where `$_` already points,
    so no path entry needs an override.
    """
    raw = yaml.safe_load((testprod / "simflow-config-l200.yaml").read_text())

    if settings_by_tier:
        shutil.copytree(
            testprod / "legend-metadata", tmp_dir / "legend-metadata", symlinks=False
        )
        for tier, overlay in settings_by_tier.items():
            f = (
                tmp_dir
                / "legend-metadata/simprod/config/tier"
                / tier
                / "l200cfg01/settings.yaml"
            )
            f.write_text(yaml.safe_dump(yaml.safe_load(f.read_text()) | overlay))
    else:
        raw["paths"]["metadata"] = str(testprod / "legend-metadata")
        raw["paths"]["config"] = str(testprod / "legend-metadata/simprod/config")

    config_path = tmp_dir / "simflow-config-l200.yaml"
    config_path.write_text(yaml.safe_dump(raw))
    return config_path


@pytest.fixture(scope="session")
def l200_config_factory():
    """Expose :func:`_l200_config` to the test modules in this directory."""
    return _l200_config


@pytest.fixture(scope="session")
def legend_currmod_paths(tmp_path_factory):
    """Run ``extract_hpge_current_pulse_model`` for both l200cfg01 runids.

    The l200cfg01 metadata has a ``default`` key so no l200data is required.
    For each runid, runs the per-detector extraction script over all modelable
    HPGes and merges the outputs into a single ``{runid}-model.yaml`` keyed by
    detector name (mirroring the ``merge_current_pulse_model_pars`` rule).
    Returns a dict mapping each runid to its merged YAML path.
    """
    out_dir = tmp_path_factory.mktemp("legend_currmod")
    config_path = _l200_config(out_dir)
    config = utils.init_simflow_context(str(config_path), workflow=None).config

    paths = {}
    for runid in _RUNIDS_L200:
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
    """Run ``make_simstat_partition_file`` with both l200cfg01 runids.

    Depends on ``legend_stp_path``; skips if remage is not installed.
    Returns the path to the output partition YAML file.
    """
    out_dir = tmp_path_factory.mktemp("legend_simstat")
    config_path = _l200_config(out_dir)
    output_file = out_dir / "partitions.yaml"

    with _override_argv(
        "make-simstat-partition-file",
        "--stp-files",
        str(legend_stp_path),
        "--runlist",
        *_RUNIDS_L200,
        "--output-file",
        str(output_file),
        "--simflow-config",
        str(config_path),
    ):
        simstat_main()

    return output_file


@pytest.fixture(scope="session")
def legend_detector_usabilities_path(tmp_path_factory):
    """Build the per-flag detector-info files for all l200cfg01 runs.

    Calls ``aggregate.gen_list_of_all_usabilities`` with the l200cfg01
    metadata (no remage or Julia required) and writes one YAML per flag
    (``usability.yaml``, ``psd_usability.yaml``,
    ``crystal_metadata_usability.yaml``). Returns the directory containing them.
    """
    out_dir = tmp_path_factory.mktemp("legend_detinfo")
    config_path = _l200_config(out_dir)

    config = utils.init_simflow_context(config_path, workflow=None).config
    detinfo = aggregate.pivot_detinfo(
        aggregate.gen_list_of_all_usabilities(config).to_dict()
    )
    for flag, mapping in detinfo.items():
        dbetto.utils.write_dict(mapping, out_dir / f"{flag}.yaml")

    return out_dir


@pytest.fixture(scope="session")
def legend_hpge_obs_paths(tmp_path_factory):
    """Run ``extract_hpge_observables_models`` for both l200cfg01 runids.

    The l200cfg01 metadata has ``default`` keys for all three observables so
    no l200data is required.  Returns a dict
    ``{runid: {eresmod: path, aoeresmod: path, psdcuts: path}}``.
    """
    out_dir = tmp_path_factory.mktemp("legend_hpge_obs")
    config_path = _l200_config(out_dir)

    paths = {}
    for runid in _RUNIDS_L200:
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
    config_path = _l200_config(out_dir)
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
        "--usability-file",
        str(legend_detector_usabilities_path / "usability.yaml"),
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
    PSD cuts, current-pulse models, and drift-time maps for the two l200cfg01
    run IDs.  Skips if remage or Julia are not installed (the stp / dtmap
    fixtures already enforce this).
    """
    out_dir = tmp_path_factory.mktemp("legend_hit")
    pars_dir = out_dir / "pars"

    for runid in _RUNIDS_L200:
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
        dtmap_dir / f"{_RUNIDS_L200[0]}-hpge-drift-time-maps.lh5",
    )
    # r001: empty dtmap → dt_map = None → PSD will be NaN
    with h5py.File(dtmap_dir / f"{_RUNIDS_L200[1]}-hpge-drift-time-maps.lh5", "w"):
        pass

    raw = yaml.safe_load((testprod / "simflow-config-l200.yaml").read_text())
    raw["paths"]["metadata"] = str(testprod / "legend-metadata")
    raw["paths"]["config"] = str(testprod / "legend-metadata/simprod/config")
    raw["paths"]["pars"] = str(pars_dir)
    config_path = out_dir / "simflow-config-l200.yaml"
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
        str(dtmap_dir / f"{_RUNIDS_L200[0]}-hpge-drift-time-maps.lh5"),
        "--currmod-files",
        str(legend_currmod_paths[_RUNIDS_L200[0]]),
        str(legend_currmod_paths[_RUNIDS_L200[1]]),
        "--simstat-part-file",
        str(legend_simstat_part_path),
        "--usability-file",
        str(legend_detector_usabilities_path / "usability.yaml"),
        "--psd-usability-file",
        str(legend_detector_usabilities_path / "psd_usability.yaml"),
        "--crystal-metadata-usability-file",
        str(legend_detector_usabilities_path / "crystal_metadata_usability.yaml"),
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
    config_path = _l200_config(out_dir)
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
        "--usability-file",
        str(legend_detector_usabilities_path / "usability.yaml"),
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
    config_path = _l200_config(out_dir)
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
