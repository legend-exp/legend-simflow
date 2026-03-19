from __future__ import annotations

from pathlib import Path

import pytest

from legendsimflow import patterns as p

SIMID = "birds_nest_K40"
JOBID = "0000"
RUNID = "l200-p02-r000-phy"
DET = "V99000A"
VOLTAGE = 4200


def test_simjob_base_segment(config):
    result = p.simjob_base_segment(config)
    assert isinstance(result, str)
    assert config.experiment in result

    result = p.simjob_base_segment(config, simid=SIMID, jobid=JOBID)
    assert SIMID in result
    assert JOBID in result


def test_log_dirname(config):
    result = p.log_dirname(config)
    assert isinstance(result, Path)
    assert config._proctime in str(result)


def test_log_filename(config):
    result = p.log_filename(config)
    assert isinstance(result, Path)

    result = p.log_filename(config, tier="stp", simid=SIMID, jobid=JOBID)
    assert isinstance(result, Path)
    assert "stp" in str(result)
    assert SIMID in str(result)
    assert JOBID in str(result)
    assert result.suffix == ".log"


def test_benchmark_filename(config):
    result = p.benchmark_filename(config)
    assert isinstance(result, Path)

    result = p.benchmark_filename(config, tier="hit", simid=SIMID, jobid=JOBID)
    assert isinstance(result, Path)
    assert "hit" in str(result)
    assert result.suffix == ".tsv"


def test_plots_dirname(config):
    for tier in ("stp", "hit", "opt", "cvt"):
        result = p.plots_dirname(config, tier)
        assert isinstance(result, Path)
        assert result.name == "plots"
        assert str(config.paths.tier[tier]) in str(result)


def test_geom_filenames(config):
    result = p.geom_config_filename(config)
    assert isinstance(result, Path)
    assert config.experiment in str(result)
    assert result.suffix == ".yaml"

    result = p.geom_gdml_filename(config)
    assert isinstance(result, Path)
    assert config.experiment in str(result)
    assert result.suffix == ".gdml"

    result = p.geom_log_filename(config)
    assert isinstance(result, Path)
    assert config.experiment in str(result)
    assert result.suffix == ".log"

    result = p.geom_gdml_filename(config, simid=SIMID, tier="stp")
    assert SIMID in str(result)
    assert "stp" in str(result)


def test_input_simjob_filename(config):
    with pytest.raises(RuntimeError, match="tier"):
        p.input_simjob_filename(config)

    # stp inputs are macro files; vtx/hit/opt inputs are lh5
    result = p.input_simjob_filename(config, tier="stp")
    assert isinstance(result, Path)
    assert result.suffix == ".mac"

    result = p.input_simjob_filename(config, tier="vtx")
    assert result.suffix == ".lh5"

    # jobid is not part of the input filename (one macro per simid)
    result = p.input_simjob_filename(config, tier="stp", simid=SIMID)
    assert isinstance(result, Path)
    assert SIMID in str(result)


def test_output_simjob_filename(config):
    with pytest.raises(RuntimeError, match="tier"):
        p.output_simjob_filename(config)

    for tier in ("vtx", "stp", "hit", "opt"):
        result = p.output_simjob_filename(config, tier=tier)
        assert isinstance(result, Path)
        assert result.suffix == ".lh5"
        assert tier in str(result)

    result = p.output_simjob_filename(config, tier="stp", simid=SIMID, jobid=JOBID)
    assert SIMID in str(result)
    assert JOBID in str(result)

    files = p.output_simjob_filename(config, tier="stp", jobid=["0000", "0001"])
    assert isinstance(files, list)
    assert len(files) == 2
    assert "0000" in files[0].name
    assert "0001" in files[1].name


def test_output_simjob_regex(config):
    with pytest.raises(RuntimeError, match="tier"):
        p.output_simjob_regex(config)

    result = p.output_simjob_regex(config, tier="stp")
    assert isinstance(result, str)
    assert config.experiment in result
    assert "stp" in result
    assert "*" in result


def test_input_simid_filenames(config):
    # input filename has no {jobid}, so n_macros does not affect the count
    result = p.input_simid_filenames(config, 3, tier="stp", simid=SIMID)
    assert isinstance(result, list)
    assert all(isinstance(f, Path) for f in result)
    assert all(SIMID in str(f) for f in result)


def test_output_simid_filenames(config):
    result = p.output_simid_filenames(config, 3, tier="stp", simid=SIMID)
    assert isinstance(result, list)
    assert len(result) == 3
    assert all(isinstance(f, Path) for f in result)
    assert all(SIMID in str(f) for f in result)


def test_vtx_filename_for_stp(config):
    result = p.vtx_filename_for_stp(config, "lar_hpge_shell_K42")
    assert isinstance(result, Path)
    assert "vtx" in str(result)

    assert p.vtx_filename_for_stp(config, "hpge_bulk_high_thr_Rn222_to_Po214") == []
    assert p.vtx_filename_for_stp(config, "phbr_surface_Ra228_to_Ac228") == []


def test_plot_tier_filenames(config):
    result = p.plot_tier_stp_vertices_filename(config)
    assert isinstance(result, Path)
    assert result.suffix == ".pdf"
    assert "stp" in str(result)
    assert result.parent.name == "plots"

    result = p.plot_tier_stp_vertices_filename(config, simid=SIMID)
    assert SIMID in result.name

    result = p.plot_tier_hit_observables_filename(config)
    assert isinstance(result, Path)
    assert "hit" in str(result)
    assert result.parent.name == "plots"

    result = p.plot_tier_opt_observables_filename(config)
    assert isinstance(result, Path)
    assert "opt" in str(result)
    assert result.parent.name == "plots"

    result = p.plot_tier_cvt_observables_filename(config)
    assert isinstance(result, Path)
    assert "cvt" in str(result)
    assert result.parent.name == "plots"


def test_dtmap_filenames(config):
    result = p.output_dtmap_filename(config, hpge_detector=DET, hpge_voltage=VOLTAGE)
    assert isinstance(result, Path)
    assert DET in str(result)
    assert f"{VOLTAGE}V" in str(result)
    assert result.suffix == ".lh5"
    assert "singles" in str(result)

    result = p.output_dtmap_merged_filename(config, runid=RUNID)
    assert isinstance(result, Path)
    assert RUNID in str(result)
    assert result.suffix == ".lh5"

    result = p.log_dtmap_filename(config, hpge_detector=DET, hpge_voltage=VOLTAGE)
    assert isinstance(result, Path)
    assert DET in str(result)
    assert result.suffix == ".log"

    result = p.plot_dtmap_filename(config, hpge_detector=DET, hpge_voltage=VOLTAGE)
    assert isinstance(result, Path)
    assert DET in str(result)
    assert f"{VOLTAGE}V" in str(result)
    assert result.suffix == ".pdf"
    assert "singles" in str(result)
    assert result.parent.name == "plots"

    result = p.benchmark_dtmap_filename(config, hpge_detector=DET, hpge_voltage=VOLTAGE)
    assert isinstance(result, Path)
    assert DET in str(result)
    assert result.suffix == ".tsv"


def test_currmod_filenames(config):
    result = p.input_currmod_evt_idx_file(config, runid=RUNID, hpge_detector=DET)
    assert isinstance(result, Path)
    assert RUNID in str(result)
    assert DET in str(result)

    result = p.output_currmod_filename(config, runid=RUNID, hpge_detector=DET)
    assert isinstance(result, Path)
    assert RUNID in str(result)
    assert DET in str(result)
    assert result.suffix == ".yaml"

    result = p.output_currmod_merged_filename(config, runid=RUNID)
    assert isinstance(result, Path)
    assert RUNID in str(result)
    assert result.suffix == ".yaml"

    result = p.log_currmod_filename(config, runid=RUNID, hpge_detector=DET)
    assert isinstance(result, Path)
    assert RUNID in str(result)
    assert DET in str(result)
    assert result.suffix == ".log"

    result = p.plot_currmod_filename(config, runid=RUNID, hpge_detector=DET)
    assert isinstance(result, Path)
    assert RUNID in str(result)
    assert DET in str(result)
    assert result.suffix == ".pdf"
    assert result.parent.name == "plots"


def test_hpge_obs_model_filenames(config):
    result = p.output_eresmod_filename(config, runid=RUNID)
    assert isinstance(result, Path)
    assert RUNID in str(result)
    assert result.suffix == ".yaml"
    assert "eresmod" in str(result)

    result = p.output_aoeresmod_filename(config, runid=RUNID)
    assert isinstance(result, Path)
    assert RUNID in str(result)
    assert result.suffix == ".yaml"
    assert "aoeresmod" in str(result)

    result = p.output_psdcuts_filename(config, runid=RUNID)
    assert isinstance(result, Path)
    assert RUNID in str(result)
    assert result.suffix == ".yaml"
    assert "psdcuts" in str(result)


def test_simstat_log_filename(config):
    result = p.log_simstat_part_filename(config, simid=SIMID)
    assert isinstance(result, Path)
    assert SIMID in str(result)
    assert result.suffix == ".log"
    assert "simstat" in str(result)


def test_cvt_filenames(config):
    result = p.tier_cvt_base_segment(config)
    assert isinstance(result, str)
    assert config.experiment in result

    result = p.tier_cvt_base_segment(config, simid=SIMID)
    assert SIMID in result

    result = p.output_tier_cvt_filename(config)
    assert isinstance(result, Path)
    assert result.suffix == ".lh5"

    result = p.output_tier_cvt_filename(config, simid=SIMID)
    assert SIMID in str(result)

    result = p.log_tier_cvt_filename(config)
    assert isinstance(result, Path)
    assert result.suffix == ".log"

    result = p.benchmark_tier_cvt_filename(config)
    assert isinstance(result, Path)
    assert result.suffix == ".tsv"


def test_pdf_filenames(config):
    result = p.pdffile_rel_basename(simid=SIMID)
    assert isinstance(result, str)
    assert SIMID in result

    result = p.output_pdf_filename(config, simid=SIMID)
    assert isinstance(result, str)
    assert SIMID in result
    assert result.endswith(".lh5")

    result = p.log_pdffile_path(config, simid=SIMID)
    assert isinstance(result, str)
    assert SIMID in result
    assert result.endswith(".log")

    result = p.benchmark_pdffile_path(config, simid=SIMID)
    assert isinstance(result, str)
    assert SIMID in result
    assert result.endswith(".tsv")
