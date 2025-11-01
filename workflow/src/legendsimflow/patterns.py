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
# TODO:
# - use coherent naming
# - use attribute access everywhere

from __future__ import annotations

from pathlib import Path

from dbetto import AttrsDict
from snakemake.io import expand

FILETYPES = AttrsDict(
    {
        "input": {
            "ver": ".mac",
            "stp": ".mac",
            "hit": ".lh5",
            "evt": ".lh5",
            "pdf": ".lh5",
        },
        "output": {
            "ver": ".lh5",
            "stp": ".lh5",
            "hit": ".lh5",
            "evt": ".lh5",
            "pdf": ".lh5",
        },
    }
)


def simjob_rel_basename(config, **kwargs):
    """Formats a partial output path for a `simid` and `jobid`."""
    return expand(
        "{simid}/" + config.experiment + "-{simid}_{jobid}",
        **kwargs,
        allow_missing=True,
    )[0]


def log_filename(config, time, **kwargs):
    """Formats a log file path for a `simid` and `jobid`."""
    pat = str(
        Path(config.paths.log)
        / time
        / "{tier}"
        / (simjob_rel_basename(config) + "-tier_{tier}.log")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


def benchmark_filename(config, **kwargs):
    """Formats a benchmark file path for a `simid` and `jobid`."""
    pat = str(
        Path(config.paths.benchmarks)
        / "{tier}"
        / (simjob_rel_basename(config) + "-tier_{tier}.tsv")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


def plots_filepath(config, **kwargs):
    """Formats a benchmark file path for a `simid` and `jobid`."""
    pat = str(Path(config.paths.plots) / "{tier}" / "{simid}")
    return expand(pat, **kwargs, allow_missing=True)[0]


# geometry


def geom_config(config, **kwargs):
    pat = str(
        Path(config.paths.geom)
        / (config.experiment + "-{simid}-tier_{tier}-geom-config.yaml")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


def geom_gdml_filename(config, **kwargs):
    pat = str(
        Path(config.paths.geom) / (config.experiment + "-{simid}-tier_{tier}-geom.gdml")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


def geom_log_filename(config, time, **kwargs):
    pat = str(
        Path(config.paths.log)
        / time
        / "geom"
        / (config.experiment + "-{simid}-tier_{tier}-geom.log")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


# ver, stp, hit tiers


def input_simjob_filename(config, **kwargs):
    """Returns the full path to the input file for a `simid`, `tier` and job index."""
    tier = kwargs.get("tier")

    if tier is None:
        msg = "the 'tier' argument is mandatory"
        raise RuntimeError(msg)

    fname = config.experiment + "-{simid}" + f"-tier_{tier}" + FILETYPES["input"][tier]
    expr = str(Path(config.paths.macros) / f"{tier}" / fname)
    return expand(expr, **kwargs, allow_missing=True)[0]


def output_simjob_filename(config, **kwargs):
    """Returns the full path to the output file for a `simid`, `tier` and job index."""
    tier = kwargs.get("tier")

    if tier is None:
        msg = "the 'tier' argument is mandatory"
        raise RuntimeError(msg)

    fname = simjob_rel_basename(config) + f"-tier_{tier}" + FILETYPES["output"][tier]
    expr = str(Path(config.paths[f"tier_{tier}"]) / fname)
    return expand(expr, **kwargs, allow_missing=True)[0]


def output_simjob_regex(config, **kwargs):
    tier = kwargs.get("tier")

    if tier is None:
        msg = "the 'tier' argument is mandatory"
        raise RuntimeError(msg)

    fname = config.experiment + "-*-tier_{tier}" + FILETYPES["output"][tier]
    expr = str(Path(config["paths"][f"tier_{tier}"]) / "{simid}" / fname)
    return expand(expr, **kwargs, allow_missing=True)[0]


def input_simid_filenames(config, n_macros, **kwargs):
    """Returns the full path to `n_macros` input files for a `simid`. Needed by
    script that generates all macros for a `simid`.
    """
    pat = input_simjob_filename(config, **kwargs)
    jobids = expand("{id:>04d}", id=list(range(n_macros)))
    return expand(pat, jobid=jobids, **kwargs, allow_missing=True)


def output_simid_filenames(config, n_macros, **kwargs):
    """Returns the full path to `n_macros` output files for a `simid`."""
    pat = output_simjob_filename(config, **kwargs)
    jobids = expand("{id:>04d}", id=list(range(n_macros)))
    return expand(pat, jobid=jobids, **kwargs, allow_missing=True)


# def ver_filename_for_stp(config, simid):
#     """Returns the vertices file needed for the 'stp' tier job, if needed. Used
#     as lambda function in the `build_tier_stp` Snakemake rule."""
#     sconfig = get_simconfig(config, "stp", simid)
#     if "vertices" in sconfig:
#         return output_simjob_filename(config, tier="ver", simid=sconfig.vertices)
#     return []


# drift time maps


def output_dtmap_filename(config, **kwargs):
    pat = str(
        Path(config.paths.dtmaps) / "{runid}-{hpge_detector}-hpge-drift-time-map.lh5"
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


def output_dtmap_merged_filename(config, **kwargs):
    pat = str(Path(config.paths.dtmaps) / "{runid}-hpge-drift-time-maps.lh5")
    return expand(pat, **kwargs, allow_missing=True)[0]


def log_dtmap_filename(config, time, **kwargs):
    pat = str(
        Path(config.paths.log)
        / time
        / "hpge_dtmaps"
        / "{runid}-{hpge_detector}-drift-time-map.log"
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


# evt tier


def evtfile_rel_basename(**kwargs):
    return expand("{simid}/{simid}_{runid}-tier_evt", **kwargs, allow_missing=True)[0]


def output_evt_filename(config, **kwargs):
    expr = str(
        Path(config["paths"]["tier_evt"])
        / (evtfile_rel_basename() + FILETYPES["output"]["evt"])
    )
    return expand(expr, **kwargs, allow_missing=True)[0]


def log_evtfile_path(config, time, **kwargs):
    pat = str(
        Path(config["paths"]["log"]) / time, "evt" / (evtfile_rel_basename() + ".log")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


def benchmark_evtfile_path(config, **kwargs):
    pat = str(
        Path(config["paths"]["benchmarks"]) / "evt" / (evtfile_rel_basename() + ".tsv")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


# pdf tier


def pdffile_rel_basename(**kwargs):
    return expand("{simid}/{simid}-tier_pdf", **kwargs, allow_missing=True)[0]


def output_pdf_filename(config, **kwargs):
    expr = str(
        Path(config["paths"]["tier_pdf"])
        / (pdffile_rel_basename() + FILETYPES["output"]["pdf"])
    )
    return expand(expr, **kwargs, allow_missing=True)[0]


def log_pdffile_path(config, time, **kwargs):
    pat = str(
        Path(config["paths"]["log"]) / time / "pdf" / (pdffile_rel_basename() + ".log")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


def benchmark_pdffile_path(config, **kwargs):
    pat = str(
        Path(config["paths"]["benchmarks"]) / "pdf" / (pdffile_rel_basename() + ".tsv")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]
