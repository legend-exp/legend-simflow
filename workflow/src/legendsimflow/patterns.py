# Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Prepare pattern strings to be used in Snakemake rules.

Extra keyword arguments are typically interpreted as variables to be
substituted in the returned (structure of) strings. They are passed to
:func:`snakemake.io.expand`.

Definitions:

- ``simid``: string identifier for the simulation run
- ``simjob``: one job of a simulation run (corresponds to one macro file and one output file)
- ``jobid``: zero-padded integer (i.e., a string) used to label a simulation job
"""

from __future__ import annotations

from pathlib import Path

from legendmeta.police import validate_dict_schema
from snakemake.io import expand

from . import SimflowConfig
from . import metadata as metautils


def _expand(pattern: str | Path, keep_list: bool = False, **kwargs) -> str | Path:
    # stringfy
    _str_pattern = pattern.as_posix() if isinstance(pattern, Path) else pattern

    exp = expand(_str_pattern, **kwargs, allow_missing=True)

    if isinstance(pattern, Path):
        exp = [Path(e) for e in exp]

    if len(exp) == 1 and not keep_list:
        exp = exp[0]

    return exp


def simjob_base_segment(config: SimflowConfig, **kwargs) -> str:
    """Formats a segment for a path including wildcards `simid` and `jobid`."""
    return _expand("{simid}/" + config.experiment + "-{simid}-job_{jobid}", **kwargs)


def log_dirname(config: SimflowConfig, time: str) -> Path:
    """Directory where log files are stored."""
    return config.paths.log / time


def log_filename(config: SimflowConfig, time: str, **kwargs) -> Path:
    """Formats a log file path for a `simid` and `jobid`."""
    pat = (
        log_dirname(config, time)
        / "{tier}"
        / (simjob_base_segment(config) + "-tier_{tier}.log")
    )
    return _expand(pat, **kwargs)


def benchmark_filename(config: SimflowConfig, **kwargs) -> Path:
    """Formats a benchmark file path for a `simid` and `jobid`."""
    pat = (
        config.paths.benchmarks
        / "{tier}"
        / (simjob_base_segment(config) + "-tier_{tier}.tsv")
    )
    return _expand(pat, **kwargs)


def plots_dirname(config: SimflowConfig, **kwargs) -> Path:
    """Formats the plots directory path for a `simid` and `tier`."""
    return _expand(config.paths.plots / "{tier}" / "{simid}", **kwargs)


# geometry


def geom_config_filename(config: SimflowConfig, **kwargs) -> Path:
    pat = config.paths.geom / (
        config.experiment + "-{simid}-tier_{tier}-geom-config.yaml"
    )
    return _expand(pat, **kwargs)


def geom_gdml_filename(config: SimflowConfig, **kwargs) -> Path:
    pat = config.paths.geom / (config.experiment + "-{simid}-tier_{tier}-geom.gdml")
    return _expand(pat, **kwargs)


def geom_log_filename(config: SimflowConfig, time: str, **kwargs) -> str:
    pat = (
        config.paths.log
        / time
        / "geom"
        / (config.experiment + "-{simid}-tier_{tier}-geom.log")
    )
    return _expand(pat, **kwargs)


# vtx, stp, opt, hit tiers


def input_simjob_filename(config: SimflowConfig, **kwargs) -> Path:
    """Returns the full path to the input file for a `simid`, `tier` and job index."""
    tier = kwargs.get("tier")

    if tier is None:
        msg = "the 'tier' argument is mandatory"
        raise RuntimeError(msg)

    ext = ".mac" if tier == "stp" else ".lh5"
    fname = config.experiment + "-{simid}" + f"-tier_{tier}" + ext
    return _expand(config.paths.macros / f"{tier}" / fname, **kwargs)


def output_simjob_filename(config: SimflowConfig, **kwargs) -> Path:
    """Returns the full path to the output file for a `simid`, `tier` and job index."""
    tier = kwargs.get("tier")

    if tier is None:
        msg = "the 'tier' argument is mandatory"
        raise RuntimeError(msg)

    fname = simjob_base_segment(config) + f"-tier_{tier}.lh5"
    return _expand(config.paths[f"tier_{tier}"] / fname, **kwargs)


def output_simjob_regex(config: SimflowConfig, **kwargs) -> str:
    tier = kwargs.get("tier")

    if tier is None:
        msg = "the 'tier' argument is mandatory"
        raise RuntimeError(msg)

    fname = config.experiment + "-*-tier_{tier}.lh5"
    expr = str(Path(config["paths"][f"tier_{tier}"]) / "{simid}" / fname)
    return _expand(expr, **kwargs)


def input_simid_filenames(config: SimflowConfig, n_macros, **kwargs) -> list[Path]:
    """Returns the full path to `n_macros` input files for a `simid`. Needed by
    script that generates all macros for a `simid`.
    """
    pat = input_simjob_filename(config, **kwargs)
    jobids = expand("{id:>04d}", id=list(range(n_macros)))
    return _expand(pat, jobid=jobids, keep_list=True, **kwargs)


def output_simid_filenames(config: SimflowConfig, n_macros, **kwargs):
    """Returns the full path to `n_macros` output files for a `simid`."""
    pat = output_simjob_filename(config, **kwargs)
    jobids = expand("{id:>04d}", id=list(range(n_macros)))
    return _expand(pat, jobid=jobids, keep_list=True, **kwargs)


def vtx_filename_for_stp(config: SimflowConfig, simid: str, **kwargs) -> Path | list:
    """Returns the vertices file needed for the 'stp' tier job, if needed.

    Used as lambda function in the `build_tier_stp` Snakemake rule.
    """
    # check the stp configuration for simid
    sconfig = metautils.get_simconfig(config, "stp", simid)

    # return vertices file only if present in the configuration
    for field in ("generator", "confinement"):
        if validate_dict_schema(
            sconfig, {field: ""}, greedy=False, verbose=False
        ) and sconfig[field].startswith("~vertices:"):
            return _expand(
                output_simjob_filename(config, tier="vtx", simid=simid), **kwargs
            )

    return []


def plot_tier_stp_vertices_filename(config: SimflowConfig, **kwargs) -> Path:
    return _expand(
        plots_dirname(config, tier="stp") / "tier-stp-vertices.pdf", **kwargs
    )


def plot_tier_hit_observables_filename(config: SimflowConfig, **kwargs) -> Path:
    return _expand(
        plots_dirname(config, tier="hit") / "tier-hit-observables.pdf", **kwargs
    )


# drift time maps


def output_dtmap_filename(config: SimflowConfig, **kwargs) -> Path:
    return _expand(
        config.paths.dtmaps / "{runid}-{hpge_detector}-hpge-drift-time-map.lh5",
        **kwargs,
    )


def output_dtmap_merged_filename(config: SimflowConfig, **kwargs) -> Path:
    return _expand(config.paths.dtmaps / "{runid}-hpge-drift-time-maps.lh5", **kwargs)


def log_dtmap_filename(config: SimflowConfig, time: str, **kwargs) -> Path:
    pat = (
        log_dirname(config, time)
        / "hpge/dtmaps/{runid}-{hpge_detector}-drift-time-map.log"
    )
    return _expand(pat, **kwargs)


def plot_dtmap_merged_filename(config: SimflowConfig, **kwargs) -> Path:
    pat = config.paths.plots / "hpge/dtmaps/{runid}-drift-time-maps.pdf"
    return _expand(pat, **kwargs)


def benchmark_dtmap_filename(config: SimflowConfig, **kwargs) -> Path:
    pat = (
        config.paths.benchmarks
        / "hpge/dtmaps/{runid}-{hpge_detector}-drift-time-map.tsv"
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


# hpge current model


def input_currmod_evt_idx_file(config: SimflowConfig, **kwargs) -> Path:
    pat = config.paths.genmeta / "hpge/currmod/{runid}-{hpge_detector}-best-evt-idx.txt"
    return expand(pat, **kwargs, allow_missing=True)[0]


def output_currmod_filename(config: SimflowConfig, **kwargs) -> Path:
    return _expand(
        config.paths.genmeta / "hpge/currmod/{runid}-{hpge_detector}-model.yaml",
        **kwargs,
    )


def output_currmod_merged_filename(config: SimflowConfig, **kwargs) -> Path:
    return _expand(
        config.paths.genmeta / "hpge/currmod/{runid}-models.yaml",
        **kwargs,
    )


def log_currmod_filename(config: SimflowConfig, time: str, **kwargs) -> Path:
    pat = log_dirname(config, time) / "hpge/currmod/{runid}-{hpge_detector}-model.log"
    return _expand(pat, **kwargs)


def plot_currmod_filename(config: SimflowConfig, **kwargs) -> Path:
    pat = config.paths.plots / "hpge/currmod/{runid}-{hpge_detector}-fit-result.pdf"
    return _expand(pat, **kwargs)


# hpge energy resolution


def output_eresmod_filename(config: SimflowConfig, **kwargs) -> Path:
    return _expand(
        config.paths.genmeta / "hpge/eresmod/{runid}-model.yaml",
        **kwargs,
    )


# hit tier


def log_simstat_part_filename(config: SimflowConfig, time: str, **kwargs) -> Path:
    pat = log_dirname(config, time) / "simstat" / "{simid}-simstat-partition.log"
    return _expand(pat, **kwargs)


# cvt tier


def tier_cvt_base_segment(config: SimflowConfig, **kwargs) -> str:
    return _expand(config.experiment + "-{simid}-tier_cvt", **kwargs)


def output_tier_cvt_filename(config: SimflowConfig, **kwargs) -> Path:
    return _expand(
        config.paths.tier_cvt / (tier_cvt_base_segment(config) + ".lh5"), **kwargs
    )


def log_tier_cvt_filename(config: SimflowConfig, time, **kwargs) -> Path:
    pat = log_dirname(config, time) / "cvt" / (tier_cvt_base_segment(config) + ".log")
    return _expand(pat, **kwargs)


def benchmark_tier_cvt_filename(config: SimflowConfig, **kwargs) -> Path:
    pat = config.paths.benchmarks / "cvt" / (tier_cvt_base_segment(config) + ".tsv")
    return _expand(pat, **kwargs)


# pdf tier


def pdffile_rel_basename(**kwargs):
    return expand("{simid}/{simid}-tier_pdf", **kwargs, allow_missing=True)[0]


def output_pdf_filename(config: SimflowConfig, **kwargs):
    expr = str(Path(config["paths"]["tier_pdf"]) / (pdffile_rel_basename() + ".lh5"))
    return expand(expr, **kwargs, allow_missing=True)[0]


def log_pdffile_path(config: SimflowConfig, time, **kwargs):
    pat = str(
        Path(config["paths"]["log"]) / time / "pdf" / (pdffile_rel_basename() + ".log")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


def benchmark_pdffile_path(config: SimflowConfig, **kwargs):
    pat = str(
        Path(config["paths"]["benchmarks"]) / "pdf" / (pdffile_rel_basename() + ".tsv")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]
