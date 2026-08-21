from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import dbetto
import lh5
import numpy as np
import pytest
from dbetto import AttrsDict
from legendmeta import LegendMetadata
from lgdo import Array, Table

from legendsimflow import hpge_pars, utils

_REPO_TEMPLATES = Path(__file__).resolve().parents[1] / "templates"


def test_get_parameter_dict():
    popt = [100, 10, 60, 0.6, 100, 0.2, 60]

    popt_dict = utils._curve_fit_popt_to_dict(popt)
    assert len(popt_dict) == 7


def test_hash_string_int():
    int_hash = utils.string_to_remage_seed("blah blah legend simflow")
    assert isinstance(int_hash, int)
    assert int_hash >= 0


def test_add_field_string():
    tab = Table(size=10)

    utils.add_field_string("string", tab, "luigi was here")

    assert isinstance(tab.string, Array)
    assert np.all(tab.string.view_as("np").astype("str") == "luigi was here")


def test_get_dataflow_config(test_l200data):
    config = utils.lookup_dataflow_config(test_l200data / "v2.1.5")

    assert isinstance(config, AttrsDict)
    assert "paths" in config
    assert isinstance(config.paths.tier, Path)
    assert "$_" not in str(config.paths.tier)


def test_check_nans_leq():
    utils.check_nans_leq([1, 2, 3, 4], "boh", 0.1, min_entries=1)
    utils.check_nans_leq([[1], [2, 3], 4], "boh", 0.1, min_entries=1)
    utils.check_nans_leq([[1], [2, 3], np.nan], "boh", 0.5, min_entries=1)

    with pytest.raises(RuntimeError):
        utils.check_nans_leq([[1], [2, 3], np.nan], "boh", 0.1, min_entries=1)

    # below min_entries: warns instead of raising
    utils.check_nans_leq([[1], [2, 3], np.nan], "boh", 0.1, min_entries=100)


def test_merge_defaults():
    """Test _merge_defaults function for recursive dictionary merging."""
    # Test basic merge
    user = {"a": 1, "b": 2}
    default = {"b": 999, "c": 3}
    result = utils._merge_defaults(user, default)
    assert result == {"a": 1, "b": 2, "c": 3}

    # Test nested merge
    user = {"outer": {"inner": 1, "keep": "user"}}
    default = {"outer": {"inner": 999, "default_key": "default"}, "top": "level"}
    result = utils._merge_defaults(user, default)
    assert result["outer"]["inner"] == 1
    assert result["outer"]["keep"] == "user"
    assert result["outer"]["default_key"] == "default"
    assert result["top"] == "level"

    # Test empty dictionaries
    assert utils._merge_defaults({}, {"a": 1}) == {"a": 1}
    assert utils._merge_defaults({"a": 1}, {}) == {"a": 1}

    # Test overwriting non-dict with dict
    user = {"key": {"nested": "value"}}
    default = {"key": "simple"}
    result = utils._merge_defaults(user, default)
    assert result["key"] == {"nested": "value"}


def test_apply_path_defaults():
    """Test apply_path_defaults fills in geom and dtmaps from pars when absent."""
    pars = Path("/prod/generated/pars")

    # both absent: should be derived from pars
    paths = {"pars": pars}
    utils.apply_path_defaults(paths)
    assert paths["geom"] == pars / "geom"
    assert paths["dtmaps"] == pars / "hpge/dtmaps"

    # both present: should be left unchanged
    custom_geom = Path("/custom/geom")
    custom_dtmaps = Path("/custom/dtmaps")
    paths = {"pars": pars, "geom": custom_geom, "dtmaps": custom_dtmaps}
    utils.apply_path_defaults(paths)
    assert paths["geom"] == custom_geom
    assert paths["dtmaps"] == custom_dtmaps

    # only one absent
    paths = {"pars": pars, "geom": custom_geom}
    utils.apply_path_defaults(paths)
    assert paths["geom"] == custom_geom
    assert paths["dtmaps"] == pars / "hpge/dtmaps"


def test_setup_logdir_link():
    """Test setup_logdir_link creates symlink to log directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        log_path = tmpdir_path / "logs"

        # Create a mock config object
        config = AttrsDict({"paths": {"log": log_path}})
        proctime = "now"
        config["_proctime"] = proctime

        # Create the target directory
        log_path.mkdir(parents=True, exist_ok=True)
        (log_path / proctime).mkdir(parents=True, exist_ok=True)

        # First call should create the symlink
        utils.setup_logdir_link(config)

        link = log_path / "latest"
        assert link.exists()
        assert link.is_symlink()
        assert link.readlink() == Path(proctime)

        # Second call should update the symlink
        new_proctime = "20240102T130000Z"
        config["_proctime"] = new_proctime
        (log_path / new_proctime).mkdir(parents=True, exist_ok=True)
        utils.setup_logdir_link(config)

        assert link.exists()
        assert link.is_symlink()
        assert link.readlink() == Path(new_proctime)


def test_init_simflow_context_returns_initialized_config(config):
    """Ensure init_simflow_context returns pre-initialized configs unchanged."""
    context = utils.init_simflow_context(config, workflow=None)
    assert context.config is config


def test_init_simflow_context_loads_from_path(tmp_path, legend_testdata):
    metadata_path = Path(legend_testdata["legend/metadata"]).resolve()
    pars_path = (tmp_path / "pars").resolve()

    raw_config = {
        "paths": {
            "metadata": str(metadata_path),
            "pars": str(pars_path),
        }
    }

    config_path = tmp_path / "simflow-config.yaml"
    dbetto.utils.write_dict(raw_config, str(config_path))

    context = utils.init_simflow_context(config_path, workflow=None)
    config = context.config
    assert isinstance(config, AttrsDict)
    assert isinstance(config.metadata, LegendMetadata)
    assert config.paths.metadata == metadata_path
    assert config.paths.pars == pars_path
    assert config.paths.geom == pars_path / "geom"
    assert config.paths.dtmaps == pars_path / "hpge/dtmaps"
    assert config._proctime


def test_hash_dict():
    """Test hash_dict produces consistent JSON representation."""
    # Test basic dictionary
    d1 = {"a": 1, "b": 2, "c": 3}
    hash1 = utils.hash_dict(d1)
    assert isinstance(hash1, str)
    assert '"a": 1' in hash1
    assert '"b": 2' in hash1

    # Test order independence
    d2 = {"c": 3, "a": 1, "b": 2}
    hash2 = utils.hash_dict(d2)
    assert hash1 == hash2

    # Test nested dictionary
    d3 = {"outer": {"inner": 1}, "list": [1, 2, 3]}
    hash3 = utils.hash_dict(d3)
    assert isinstance(hash3, str)
    assert "outer" in hash3

    # Test AttrsDict
    attrs_dict = AttrsDict({"key": "value", "num": 42})
    hash4 = utils.hash_dict(attrs_dict)
    assert isinstance(hash4, str)
    assert "key" in hash4


@pytest.fixture
def tier_test_data(tmp_path):
    """Create a temporary test data directory with tier configuration."""
    # Create config YAML with setups.l200 structure
    config_file = tmp_path / "config.yaml"
    config_content = {
        "setups": {
            "l200": {
                "paths": {
                    "tier_hit": "$_/generated/tier/hit",
                    "tier_pht": "$_/generated/tier/pht",
                    "tier_raw": "$_/generated/tier/raw",
                    "par": "$_/generated/par",
                    "par_hit": "$_/generated/par/hit",
                    "par_pht": "$_/generated/par/pht",
                }
            }
        }
    }
    dbetto.utils.write_dict(config_content, str(config_file))

    # Create tier directories
    hit_dir = tmp_path / "generated" / "tier" / "hit" / "phy" / "p00" / "r000"
    raw_dir = tmp_path / "generated" / "tier" / "raw" / "phy" / "p00" / "r000"
    pht_dir = tmp_path / "generated" / "tier" / "pht" / "phy" / "p00" / "r000"

    hit_dir.mkdir(parents=True, exist_ok=True)
    pht_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    # write a few files
    lh5.write(Table(), "hit", hit_dir / "file1.lh5")
    lh5.write(Table(), "hit", hit_dir / "file2.lh5")

    lh5.write(Table(), "raw", raw_dir / "file1.lh5")
    lh5.write(Table(), "raw", raw_dir / "file2.lh5")

    # Create par directories
    par_dir = tmp_path / "generated" / "par"
    par_hit_dir = tmp_path / "generated" / "par" / "hit"
    par_pht_dir = tmp_path / "generated" / "par" / "pht"
    par_dir.mkdir(parents=True, exist_ok=True)
    par_hit_dir.mkdir(parents=True, exist_ok=True)
    par_pht_dir.mkdir(parents=True, exist_ok=True)

    return tmp_path


@pytest.fixture
def tier_test_data_direct(tmp_path):
    """Create a temporary test data directory with direct paths format (no setups.l200)."""
    # Create config YAML with direct paths structure (no setups wrapper)
    config_file = tmp_path / "config.yaml"
    config_content = {
        "paths": {
            "tier_hit": "$_/generated/tier/hit",
            "tier_pht": "$_/generated/tier/pht",
            "tier_raw": "$_/generated/tier/raw",
            "par": "$_/generated/par",
            "par_hit": "$_/generated/par/hit",
            "par_pht": "$_/generated/par/pht",
        }
    }
    dbetto.utils.write_dict(config_content, str(config_file))

    # Create tier directories
    hit_dir = tmp_path / "generated" / "tier" / "hit" / "phy" / "p00" / "r000"
    raw_dir = tmp_path / "generated" / "tier" / "raw" / "phy" / "p00" / "r000"
    pht_dir = tmp_path / "generated" / "tier" / "pht" / "phy" / "p00" / "r000"

    hit_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    pht_dir.mkdir(parents=True, exist_ok=True)

    lh5.write(Table(), "hit", hit_dir / "file1.lh5")
    lh5.write(Table(), "hit", hit_dir / "file2.lh5")

    lh5.write(Table(), "raw", raw_dir / "file1.lh5")
    lh5.write(Table(), "raw", raw_dir / "file2.lh5")

    # Create par directories
    par_dir = tmp_path / "generated" / "par"
    par_hit_dir = tmp_path / "generated" / "par" / "hit"
    par_pht_dir = tmp_path / "generated" / "par" / "pht"

    par_dir.mkdir(parents=True, exist_ok=True)
    par_hit_dir.mkdir(parents=True, exist_ok=True)
    par_pht_dir.mkdir(parents=True, exist_ok=True)

    return tmp_path


def test_lookup_file_paths(tier_test_data, tier_test_data_direct):
    for path in [tier_test_data, tier_test_data_direct]:
        for dtype in ["raw", "hit"]:
            assert set(
                hpge_pars.lookup_file_paths(
                    path, "l200-p00-r000-phy", hit_tier_name="hit"
                )[dtype]
            ) == {
                path
                / "generated"
                / "tier"
                / dtype
                / "phy"
                / "p00"
                / "r000"
                / "file1.lh5",
                path
                / "generated"
                / "tier"
                / dtype
                / "phy"
                / "p00"
                / "r000"
                / "file2.lh5",
            }


def test_lookup_dataflow_config_with_fixture(tier_test_data):
    """Test lookup_dataflow_config with a proper fixture."""
    config = utils.lookup_dataflow_config(tier_test_data)

    assert isinstance(config, AttrsDict)
    assert "paths" in config
    assert "tier_hit" in config.paths
    assert "tier_pht" in config.paths
    # Verify variable substitution happened
    assert "$_" not in str(config.paths.tier_hit)
    assert "generated/tier/hit" in str(config.paths.tier_hit)


def test_get_hit_tier_name_with_pht(tier_test_data):
    """Test get_hit_tier_name returns 'pht' when both tiers exist."""
    # When both hit and pht exist, pht should be preferred
    tier_name = utils.get_hit_tier_name(str(tier_test_data))
    assert tier_name == "pht"


def test_get_hit_tier_name_with_hit_only(tmp_path):
    """Test get_hit_tier_name returns 'hit' when only hit tier exists."""
    # Create config with only hit tier
    config_content = """
paths:
  tier_hit: $_/generated/tier/hit
"""
    config_file = tmp_path / "config.yaml"
    config_file.write_text(config_content)

    # Create only hit directory
    hit_dir = tmp_path / "generated" / "tier" / "hit"
    hit_dir.mkdir(parents=True, exist_ok=True)

    tier_name = utils.get_hit_tier_name(str(tmp_path))
    assert tier_name == "hit"


def test_get_hit_tier_name_missing_dirs(tmp_path):
    """Test get_hit_tier_name raises error when directories don't exist."""
    # Create config but no actual directories
    config_content = """
paths:
  tier_hit: $_/generated/tier/hit
  tier_pht: $_/generated/tier/pht
"""
    config_file = tmp_path / "config.yaml"
    config_file.write_text(config_content)

    with pytest.raises(RuntimeError, match="does not contain a valid pht or hit tier"):
        utils.get_hit_tier_name(str(tmp_path))


def test_init_generated_pars_db(tier_test_data):
    """Test init_generated_pars_db initializes pars database correctly."""
    # Test getting the full par database
    par_db = utils.init_generated_pars_db(tier_test_data, tier=None, lazy=True)
    assert par_db is not None
    # Verify the path is correct (TextDB has __path__ attribute)
    assert "generated/par" in repr(par_db)

    # Test getting the hit tier pars database
    par_hit_db = utils.init_generated_pars_db(tier_test_data, tier="hit", lazy=True)
    assert par_hit_db is not None
    assert "generated/par/hit" in repr(par_hit_db)

    # Test getting the pht tier pars database
    par_pht_db = utils.init_generated_pars_db(tier_test_data, tier="pht", lazy=True)
    assert par_pht_db is not None
    assert "generated/par/pht" in repr(par_pht_db)


def test_lookup_dataflow_config_direct_format(tier_test_data_direct):
    """Test lookup_dataflow_config with direct paths format (no setups.l200)."""
    config = utils.lookup_dataflow_config(tier_test_data_direct)

    assert isinstance(config, AttrsDict)
    assert "paths" in config
    assert "tier_hit" in config.paths
    assert "tier_pht" in config.paths
    # Verify variable substitution happened
    assert "$_" not in str(config.paths.tier_hit)
    assert "generated/tier/hit" in str(config.paths.tier_hit)


def test_get_hit_tier_name_direct_format(tier_test_data_direct):
    """Test get_hit_tier_name with direct paths format."""
    tier_name = utils.get_hit_tier_name(str(tier_test_data_direct))
    assert tier_name == "pht"


def test_init_generated_pars_db_direct_format(tier_test_data_direct):
    """Test init_generated_pars_db with direct paths format."""
    # Test getting the full par database
    par_db = utils.init_generated_pars_db(tier_test_data_direct, tier=None, lazy=True)
    assert par_db is not None
    assert "generated/par" in repr(par_db)

    # Test getting the hit tier pars database
    par_hit_db = utils.init_generated_pars_db(
        tier_test_data_direct, tier="hit", lazy=True
    )
    assert par_hit_db is not None
    assert "generated/par/hit" in repr(par_hit_db)

    # Test getting the pht tier pars database
    par_pht_db = utils.init_generated_pars_db(
        tier_test_data_direct, tier="pht", lazy=True
    )
    assert par_pht_db is not None
    assert "generated/par/pht" in repr(par_pht_db)


def test_sanitize_dict():
    defaults = {"a": {"x": -1.5, "y": 3}, "b": {"flag": True, "name": "ok"}}
    read = {"a": {"x": "nope", "y": 4.2}, "b": {"flag": "yes"}}

    out = utils.sanitize_dict_with_defaults(read, defaults)

    assert out["a"]["x"] == -1.5  # replaced (non-numeric)
    assert out["a"]["y"] == 4.2  # kept (numeric)
    assert out["b"]["flag"] is True  # replaced (wrong type)
    assert out["b"]["name"] == "ok"  # filled from default


def test_get_dict_value():
    d = {"a": {"b": 3, "c": 2}}
    assert utils.get_dict_value(d, "a") == {"b": 3, "c": 2}
    assert utils.get_dict_value(d, "a.b") == 3
    assert utils.get_dict_value(d, "a.b.d") is None
    assert utils.get_dict_value(d, "a.b.d", 69) == 69


@pytest.fixture
def link_setup(tmp_path, monkeypatch):
    """Fake simflow source repo + a prod cycle as cwd + a config of all-default paths."""
    repo = tmp_path / "simflow"
    workflow_basedir = repo / "workflow"
    workflow_basedir.mkdir(parents=True)
    (repo / "templates").mkdir()
    shutil.copy(_REPO_TEMPLATES / "default.yaml", repo / "templates" / "default.yaml")

    cycle = tmp_path / "cycle"
    cycle.mkdir()
    monkeypatch.chdir(cycle)

    config = AttrsDict(
        {
            "paths": AttrsDict(
                {
                    "pars": cycle / "generated/pars",
                    "macros": cycle / "generated/macros",
                    "geom": cycle / "generated/pars/geom",
                    "dtmaps": cycle / "generated/pars/hpge/dtmaps",
                    "tier": AttrsDict(
                        {
                            t: cycle / f"generated/tier/{t}"
                            for t in (
                                "vtx",
                                "stp",
                                "par",
                                "opt",
                                "hit",
                                "evt",
                                "cvt",
                                "pdf",
                            )
                        }
                    ),
                }
            )
        }
    )
    return workflow_basedir, cycle, config


def test_link_external_paths_links_overrides(link_setup, tmp_path):
    workflow_basedir, cycle, config = link_setup

    ext_hit = tmp_path / "external" / "hit"
    ext_opt = tmp_path / "external" / "opt"
    ext_pars = tmp_path / "external" / "pars"
    for d in (ext_hit, ext_opt, ext_pars):
        d.mkdir(parents=True)
    config.paths.tier["hit"] = ext_hit
    config.paths.tier["opt"] = ext_opt
    config.paths["pars"] = ext_pars

    utils.link_external_paths(config, workflow_basedir)

    for default, target in (
        (cycle / "generated/tier/hit", ext_hit),
        (cycle / "generated/tier/opt", ext_opt),
        (cycle / "generated/pars", ext_pars),
    ):
        assert default.is_symlink()
        assert not default.readlink().is_absolute()
        assert default.resolve() == target.resolve()

    assert not (cycle / "generated/tier/evt").exists()


def test_link_external_paths_noop_when_at_default(link_setup):
    workflow_basedir, cycle, config = link_setup

    utils.link_external_paths(config, workflow_basedir)

    for tier in ("hit", "opt", "evt"):
        assert not (cycle / f"generated/tier/{tier}").exists()
    assert not (cycle / "generated/pars").exists()


def test_link_external_paths_skips_real_default_directories(link_setup, tmp_path):
    workflow_basedir, cycle, config = link_setup

    ext_hit = tmp_path / "external" / "hit"
    ext_hit.mkdir(parents=True)
    config.paths.tier["hit"] = ext_hit

    real = cycle / "generated/tier/hit"
    real.mkdir(parents=True)
    (real / "data.lh5").touch()

    utils.link_external_paths(config, workflow_basedir)

    assert not real.is_symlink()
    assert (real / "data.lh5").is_file()


def test_link_external_paths_is_idempotent(link_setup, tmp_path):
    workflow_basedir, cycle, config = link_setup

    ext_hit = tmp_path / "external" / "hit"
    ext_hit.mkdir(parents=True)
    config.paths.tier["hit"] = ext_hit

    utils.link_external_paths(config, workflow_basedir)
    default = cycle / "generated/tier/hit"
    target1 = default.readlink()

    utils.link_external_paths(config, workflow_basedir)
    target2 = default.readlink()

    assert target1 == target2


def test_link_external_paths_replaces_wrong_symlink(link_setup, tmp_path):
    workflow_basedir, cycle, config = link_setup

    ext_hit = tmp_path / "external" / "hit"
    ext_hit.mkdir(parents=True)
    config.paths.tier["hit"] = ext_hit

    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    default = cycle / "generated/tier/hit"
    default.parent.mkdir(parents=True)
    default.symlink_to(elsewhere)

    utils.link_external_paths(config, workflow_basedir)

    assert default.is_symlink()
    assert default.resolve() == ext_hit.resolve()


def test_link_external_paths_removes_stale_symlink(link_setup, tmp_path):
    workflow_basedir, cycle, config = link_setup

    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    default = cycle / "generated/tier/hit"
    default.parent.mkdir(parents=True)
    default.symlink_to(elsewhere)

    utils.link_external_paths(config, workflow_basedir)

    assert not default.is_symlink()
    assert not default.exists()


def test_link_external_paths_geom_dtmaps_default(link_setup, tmp_path):
    workflow_basedir, cycle, config = link_setup
    del config.paths["pars"]
    ext_geom = tmp_path / "external" / "geom"
    ext_dtmaps = tmp_path / "external" / "dtmaps"
    ext_geom.mkdir(parents=True)
    ext_dtmaps.mkdir(parents=True)
    config.paths["geom"] = ext_geom
    config.paths["dtmaps"] = ext_dtmaps

    utils.link_external_paths(config, workflow_basedir)

    default_geom = cycle / "generated/pars/geom"
    default_dtmaps = cycle / "generated/pars/hpge/dtmaps"
    assert default_geom.is_symlink()
    assert default_geom.resolve() == ext_geom.resolve()
    assert default_dtmaps.is_symlink()
    assert default_dtmaps.resolve() == ext_dtmaps.resolve()


def test_link_external_paths_strips_dvs_ro_prefix(link_setup, tmp_path):
    """Strip the NERSC /dvs_ro mount prefix before symlinking.

    Override paths under the NERSC /dvs_ro/... read-only mount must be
    normalized to /global/... before computing the relative symlink — otherwise
    the link traverses all the way to / and back down through the mount.
    """
    workflow_basedir, cycle, config = link_setup

    ext_hit = tmp_path / "external" / "hit"
    ext_hit.mkdir(parents=True)
    # the user's config points at the /dvs_ro mirror of the canonical /global path
    config.paths.tier["hit"] = Path("/dvs_ro") / ext_hit.relative_to("/")

    utils.link_external_paths(config, workflow_basedir)

    default = cycle / "generated/tier/hit"
    assert default.is_symlink()
    target = default.readlink()
    assert not target.is_absolute()
    assert "dvs_ro" not in str(target)


def test_link_external_paths_creates_intermediate_dirs(link_setup, tmp_path):
    workflow_basedir, cycle, config = link_setup

    ext_hit = tmp_path / "external" / "hit"
    ext_hit.mkdir(parents=True)
    config.paths.tier["hit"] = ext_hit

    assert not (cycle / "generated").exists()

    utils.link_external_paths(config, workflow_basedir)

    assert (cycle / "generated/tier/hit").is_symlink()
