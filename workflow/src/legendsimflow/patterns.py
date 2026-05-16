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
    """Expand a path pattern with Snakemake wildcards.

    Returning a scalar unless `keep_list` is set.
    """
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


def log_dirname(config: SimflowConfig) -> Path:
    """Directory where log files are stored."""
    return config.paths.log / config._proctime


def log_filename(config: SimflowConfig, **kwargs) -> Path:
    """Formats a log file path for a `simid` and `jobid`."""
    pat = (
        log_dirname(config)
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


def plots_dirname(config: SimflowConfig, tier: str) -> Path:
    """Returns the plots directory path for a `tier`."""
    return config.paths.tier[tier] / "plots"


def plots_tarball_filename(config: SimflowConfig) -> Path:
    """The path to the plots archive tarball for the current production cycle.

    The Simflow has no explicit knowledge of the production cycle name, so the
    name of the directory where the Simflow lives is used as a proxy.
    """
    cycle_dir = config.paths.generated.parent
    return config.paths.generated / "tarballs" / (cycle_dir.name + "-plots.tar.xz")


def pdf_tarball_filename(config: SimflowConfig) -> Path:
    """The path to the pdf tier archive tarball for the current production cycle.

    The Simflow has no explicit knowledge of the production cycle name, so the
    name of the directory where the Simflow lives is used as a proxy.
    """
    cycle_dir = config.paths.generated.parent
    return config.paths.generated / "tarballs" / (cycle_dir.name + "-pdfs.tar.xz")


# geometry


def geom_config_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the geometry configuration YAML file for a `tier` and `simid`."""
    pat = config.paths.geom / (
        config.experiment + "-{simid}-tier_{tier}-geom-config.yaml"
    )
    return _expand(pat, **kwargs)


def geom_gdml_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the GDML geometry file for a `tier` and `simid`."""
    pat = config.paths.geom / (config.experiment + "-{simid}-tier_{tier}-geom.gdml")
    return _expand(pat, **kwargs)


def geom_log_filename(config: SimflowConfig, **kwargs) -> str:
    """The log file path for geometry generation for a `tier` and `simid`."""
    pat = (
        config.paths.log
        / config._proctime
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
    return _expand(config.paths.macros / fname, **kwargs)


def output_simjob_filename(config: SimflowConfig, **kwargs) -> Path:
    """Returns the full path to the output file for a `simid`, `tier` and job index."""
    tier = kwargs.get("tier")

    if tier is None:
        msg = "the 'tier' argument is mandatory"
        raise RuntimeError(msg)

    fname = simjob_base_segment(config) + f"-tier_{tier}.lh5"
    return _expand(config.paths.tier[tier] / fname, **kwargs)


def output_simjob_regex(config: SimflowConfig, **kwargs) -> str:
    """A glob-style regex matching all output files for a `tier`."""
    tier = kwargs.get("tier")

    if tier is None:
        msg = "the 'tier' argument is mandatory"
        raise RuntimeError(msg)

    fname = config.experiment + "-*-tier_{tier}.lh5"
    expr = str(config.paths.tier[tier] / "{simid}" / fname)
    return _expand(expr, **kwargs)


def input_simid_filenames(config: SimflowConfig, n_macros, **kwargs) -> list[Path]:
    """Returns the full path to `n_macros` input files for a `simid`.

    Needed by script that generates all macros for a `simid`.
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
    """The path to the primary vertex validation plot for a `stp` `simid`."""
    return _expand(
        plots_dirname(config, "stp") / "{simid}-tier-stp-vertices.pdf", **kwargs
    )


def plot_tier_hit_observables_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the observable validation plot for a `hit` `simid`."""
    return _expand(
        plots_dirname(config, "hit") / "{simid}-tier-hit-observables.pdf", **kwargs
    )


def plot_tier_opt_observables_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the observable validation plot for an `opt` `simid`."""
    return _expand(
        plots_dirname(config, "opt") / "{simid}-tier-opt-observables.pdf", **kwargs
    )


def plot_tier_cvt_observables_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the observable validation plot for a `cvt` `simid`."""
    return _expand(
        plots_dirname(config, "cvt") / "{simid}-tier-cvt-observables.pdf", **kwargs
    )


# drift time maps


def output_dtmap_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the HPGe drift time map file for a detector and voltage."""
    return _expand(
        config.paths.dtmaps
        / "singles/{hpge_detector}-{hpge_voltage}V-hpge-drift-time-map.lh5",
        **kwargs,
    )


def output_dtmap_merged_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the merged HPGe drift time map file for a `runid`."""
    return _expand(config.paths.dtmaps / "{runid}-hpge-drift-time-maps.lh5", **kwargs)


def output_psl_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the HPGe psl file for a detector and voltage."""
    return _expand(
        config.paths.pulse_shape_lib
        / "singles/{hpge_detector}-{hpge_voltage}V-hpge-pulse-shape-library.lh5",
        **kwargs,
    )


def output_psl_merged_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the merged HPGe psl file for a `runid`."""
    return _expand(
        config.paths.pulse_shape_lib / "{runid}-hpge-pulse-shape-library.lh5", **kwargs
    )


def log_psl_filename(config: SimflowConfig, **kwargs) -> Path:
    """The log file path for drift time map generation for a detector and voltage."""
    pat = (
        log_dirname(config)
        / "hpge/pulse_shape_lib/{hpge_detector}-{hpge_voltage}V-drift-time-map.log"
    )
    return _expand(pat, **kwargs)


def plot_psl_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the pulse shape library validation plot for a detector and voltage."""
    pat = (
        config.paths.pulse_shape_lib
        / "singles/plots/{hpge_detector}-{hpge_voltage}V-drift-time-map.pdf"
    )
    return _expand(pat, **kwargs)


def benchmark_psl_filename(config: SimflowConfig, **kwargs) -> Path:
    """The benchmark file path for pulse shape library generation for a detector and voltage."""
    pat = (
        config.paths.benchmarks
        / "hpge/pulse_shape_lib/{hpge_detector}-{hpge_voltage}V-drift-time-map.tsv"
    )
    return _expand(pat, **kwargs)


def log_dtmap_filename(config: SimflowConfig, **kwargs) -> Path:
    """The log file path for drift time map generation for a detector and voltage."""
    pat = (
        log_dirname(config)
        / "hpge/dtmaps/{hpge_detector}-{hpge_voltage}V-drift-time-map.log"
    )
    return _expand(pat, **kwargs)


def plot_dtmap_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the drift time map validation plot for a detector and voltage."""
    pat = (
        config.paths.dtmaps
        / "singles/plots/{hpge_detector}-{hpge_voltage}V-drift-time-map.pdf"
    )
    return _expand(pat, **kwargs)


def benchmark_dtmap_filename(config: SimflowConfig, **kwargs) -> Path:
    """The benchmark file path for drift time map generation for a detector and voltage."""
    pat = (
        config.paths.benchmarks
        / "hpge/dtmaps/{hpge_detector}-{hpge_voltage}V-drift-time-map.tsv"
    )
    return _expand(pat, **kwargs)


# hpge pulse shape libraries


def output_ideal_psl_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the ideal HPGe pulse shape library for a detector and voltage."""
    pat = (
        config.paths.pars
        / "hpge/psl/ideal/singles/{hpge_detector}-{hpge_voltage}V-ideal-pulse-shape-lib.lh5"
    )
    return _expand(pat, **kwargs)


def log_ideal_psl_filename(config: SimflowConfig, **kwargs) -> Path:
    """The log file path for ideal pulse shape library generation for a detector and voltage."""
    pat = (
        log_dirname(config)
        / "hpge/psl/ideal/{hpge_detector}-{hpge_voltage}V-ideal-pulse-shape-lib.log"
    )
    return _expand(pat, **kwargs)


def benchmark_ideal_psl_filename(config: SimflowConfig, **kwargs) -> Path:
    """The benchmark file path for ideal pulse shape library generation for a detector and voltage."""
    pat = (
        config.paths.benchmarks
        / "hpge/psl/ideal/{hpge_detector}-{hpge_voltage}V-ideal-pulse-shape-lib.tsv"
    )
    return _expand(pat, **kwargs)


def output_realistic_psl_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the realistic HPGe pulse shape library for a detector and run."""
    pat = (
        config.paths.pars
        / "hpge/psl/realistic/{runid}-{hpge_detector}-pulse-shape-lib.lh5"
    )
    return _expand(pat, **kwargs)


def log_realistic_psl_filename(config: SimflowConfig, **kwargs) -> Path:
    """The log file path for realistic pulse shape library generation for a detector and run."""
    pat = (
        log_dirname(config)
        / "hpge/psl/realistic/{runid}-{hpge_detector}-pulse-shape-lib.log"
    )
    return _expand(pat, **kwargs)


def benchmark_realistic_psl_filename(config: SimflowConfig, **kwargs) -> Path:
    """The benchmark file path for realistic pulse shape library generation for a detector and run."""
    pat = (
        config.paths.benchmarks
        / "hpge/psl/realistic/{runid}-{hpge_detector}-pulse-shape-lib.tsv"
    )
    return _expand(pat, **kwargs)


# hpge current model


def input_currmod_evt_idx_file(config: SimflowConfig, **kwargs) -> Path:
    """The path to the event index file used to extract current pulse waveforms."""
    pat = config.paths.pars / "hpge/currmod/{runid}-{hpge_detector}-best-evt-idx.txt"
    return _expand(pat, **kwargs)


def output_currmod_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the per-detector HPGe current pulse model parameter file."""
    return _expand(
        config.paths.pars / "hpge/currmod/{runid}-{hpge_detector}-model.yaml",
        **kwargs,
    )


def output_currmod_merged_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the merged HPGe current pulse model parameter file for a `runid`."""
    return _expand(
        config.paths.pars / "hpge/currmod/{runid}-model.yaml",
        **kwargs,
    )


def log_currmod_filename(config: SimflowConfig, **kwargs) -> Path:
    """The log file path for current pulse model extraction for a detector and `runid`."""
    pat = log_dirname(config) / "hpge/currmod/{runid}-{hpge_detector}-model.log"
    return _expand(pat, **kwargs)


def plot_currmod_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the current pulse model fit validation plot for a detector and `runid`."""
    pat = (
        config.paths.pars / "hpge/currmod/plots/{runid}-{hpge_detector}-fit-result.pdf"
    )
    return _expand(pat, **kwargs)


# hpge energy resolution


def log_eresmod_filename(config: SimflowConfig, **kwargs) -> Path:
    """The log file path for HPGe observables model extraction for a `runid`."""
    pat = log_dirname(config) / "hpge/eresmod/{runid}-model.log"
    return _expand(pat, **kwargs)


def output_eresmod_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the HPGe energy resolution model parameter file for a `runid`."""
    return _expand(
        config.paths.pars / "hpge/eresmod/{runid}-model.yaml",
        **kwargs,
    )


def output_aoeresmod_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the HPGe A/E resolution model parameter file for a `runid`."""
    return _expand(
        config.paths.pars / "hpge/aoeresmod/{runid}-model.yaml",
        **kwargs,
    )


def output_psdcuts_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the HPGe PSD cut values file for a `runid`."""
    return _expand(
        config.paths.pars / "hpge/psdcuts/{runid}-psd-cuts.yaml",
        **kwargs,
    )


# hit tier


def simstat_part_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the simulation event statistics partitioning file."""
    pat = config.paths.pars / "simstat" / "partitions_{simid}.yaml"
    return _expand(pat, **kwargs)


def log_simstat_part_filename(config: SimflowConfig, **kwargs) -> Path:
    """The log file path for simulation event statistics partitioning for a `simid`."""
    pat = log_dirname(config) / "simstat" / "{simid}-simstat-partition.log"
    return _expand(pat, **kwargs)


# cvt tier


def tier_cvt_base_segment(config: SimflowConfig, **kwargs) -> str:
    """The base filename segment for `cvt` tier files for a `simid`."""
    return _expand(config.experiment + "-{simid}-tier_cvt", **kwargs)


def output_tier_cvt_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the merged `cvt` tier output file for a `simid`."""
    return _expand(
        config.paths.tier.cvt / (tier_cvt_base_segment(config) + ".lh5"), **kwargs
    )


def log_tier_cvt_filename(config: SimflowConfig, **kwargs) -> Path:
    """The log file path for the `cvt` tier build for a `simid`."""
    pat = log_dirname(config) / "cvt" / (tier_cvt_base_segment(config) + ".log")
    return _expand(pat, **kwargs)


def benchmark_tier_cvt_filename(config: SimflowConfig, **kwargs) -> Path:
    """The benchmark file path for the `cvt` tier build for a `simid`."""
    pat = config.paths.benchmarks / "cvt" / (tier_cvt_base_segment(config) + ".tsv")
    return _expand(pat, **kwargs)


# pdf tier


def tier_pdf_base_segment(config: SimflowConfig, **kwargs) -> str:
    """The base filename segment for `pdf` tier files for a `simid`."""
    return _expand(config.experiment + "-{simid}-tier_pdf", **kwargs)


def output_tier_pdf_filename(config: SimflowConfig, **kwargs) -> Path:
    """The path to the merged `pdf` tier output file for a `simid`."""
    return _expand(
        config.paths.tier.pdf / (tier_pdf_base_segment(config) + ".lh5"), **kwargs
    )


def log_tier_pdf_filename(config: SimflowConfig, **kwargs) -> Path:
    """The log file path for the `pdf` tier build for a `simid`."""
    pat = log_dirname(config) / "pdf" / (tier_pdf_base_segment(config) + ".log")
    return _expand(pat, **kwargs)


def benchmark_tier_pdf_filename(config: SimflowConfig, **kwargs) -> Path:
    """The benchmark file path for the `pdf` tier build for a `simid`."""
    pat = config.paths.benchmarks / "pdf" / (tier_pdf_base_segment(config) + ".tsv")
    return _expand(pat, **kwargs)
