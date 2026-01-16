from __future__ import annotations

from pathlib import Path

from legendsimflow import patterns as p


def test_all(config):
    assert isinstance(p.simjob_base_segment(config), str)
    assert isinstance(p.log_filename(config, "now"), Path)
    assert isinstance(p.plots_dirname(config), Path)
    assert isinstance(p.benchmark_filename(config), Path)
    assert isinstance(p.geom_gdml_filename(config), Path)
    assert isinstance(p.input_simjob_filename(config, tier="stp"), Path)
    assert isinstance(p.output_simjob_filename(config, tier="stp"), Path)
    assert isinstance(p.output_simjob_regex(config, tier="stp"), str)
    assert isinstance(p.input_simid_filenames(config, 2, tier="stp"), list)
    assert all(
        isinstance(p, Path) for p in p.input_simid_filenames(config, 2, tier="stp")
    )
    assert isinstance(p.output_simid_filenames(config, 2, tier="stp"), list)

    assert isinstance(p.output_dtmap_filename(config, hpge_detector="boh"), Path)
    assert isinstance(p.output_currmod_filename(config, hpge_detector="boh"), Path)

    assert isinstance(p.vtx_filename_for_stp(config, "lar_hpge_shell_K42"), Path)
    assert p.vtx_filename_for_stp(config, "hpge_bulk_high_thr_Rn222_to_Po214") == []
    assert p.vtx_filename_for_stp(config, "phbr_surface_Ra228_to_Ac228") == []

    assert isinstance(p.tier_cvt_base_segment(config), str)
    assert isinstance(p.output_tier_cvt_filename(config), Path)
    assert isinstance(p.log_tier_cvt_filename(config, "now"), Path)
