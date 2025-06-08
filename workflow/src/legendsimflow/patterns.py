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
- ``simjob``: one job of a simulation run (corresponds to one macro file and one
  output file)
- ``jobid``: zero-padded integer (i.e., a string) used to label a simulation job
"""
from __future__ import annotations

from pathlib import Path

import yaml
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


def simjob_rel_basename(**kwargs):
    """Formats a partial output path for a `simid` and `jobid`."""
    return expand("{simid}/{simid}_{jobid}", **kwargs, allow_missing=True)[0]


def run_command(config, tier):
    """Returns command to build files in tier `tier` prefixed by environment."""
    return " ".join(config["execenv"]) + " " + config["runcmd"][tier]


def log_file_path(config, time, **kwargs):
    """Formats a log file path for a `simid` and `jobid`."""
    pat = str(
        Path(config["paths"]["log"])
        / time
        / "{tier}"
        / (simjob_rel_basename() + "-tier_{tier}.log")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


def benchmark_file_path(config, **kwargs):
    """Formats a benchmark file path for a `simid` and `jobid`."""
    pat = str(
        Path(config["paths"]["benchmarks"])
        / "{tier}"
        / (simjob_rel_basename() + "-tier_{tier}.tsv")
    )
    return expand(pat, **kwargs, allow_missing=True)[0]


def plots_file_path(config, **kwargs):
    """Formats a benchmark file path for a `simid` and `jobid`."""
    pat = str(Path(config["paths"]["plots"]) / "{tier}" / "{simid}")
    return expand(pat, **kwargs, allow_missing=True)[0]


def genmacro_log_file_path(config, time, **kwargs):
    """Formats a log file path for a `simid` and `jobid`."""
    return expand(
        str(
            Path(config["paths"]["log"])
            / time
            / "macros"
            / "{tier}"
            / (simjob_rel_basename() + "-tier_{tier}.log")
        ),
        **kwargs,
        allow_missing=True,
    )[0]


def template_macro_dir(config, **kwargs):
    """Returns the directory path to the macro templates for the current `tier`."""
    tier = expand("{tier}", **kwargs, allow_missing=True)[0]
    return Path(config["paths"]["config"]) / "tier" / tier / config["experiment"]


# ver, stp, hit tiers


def macro_gen_inputs(config, tier, simid, **kwargs):
    """Return inputs for the Snakemake rules that generate macros."""
    tdir = template_macro_dir(config, tier=tier)

    with (tdir / "simconfig.yaml").open() as f:
        sconfig = yaml.safe_load(f)[simid]

    if "template" not in sconfig:
        msg = "simconfig.yaml blocks must define a 'template' field."
        raise RuntimeError(msg)

    expr = {
        "template": str(tdir / sconfig["template"]),
        "cfgfile": str(tdir / "simconfig.yaml"),
    }
    for k, v in expr.items():
        expr[k] = expand(v, **kwargs, allow_missing=True)[0]
    return expr


def input_simjob_filename(config, **kwargs):
    """Returns the full path to the input file for a `simid`, `tier` and job index."""
    tier = kwargs.get("tier")

    if tier is None:
        msg = "the 'tier' argument is mandatory"
        raise RuntimeError(msg)

    fname = simjob_rel_basename() + f"-tier_{tier}" + FILETYPES["input"][tier]
    expr = str(Path(config["paths"]["macros"]) / f"{tier}" / fname)
    return expand(expr, **kwargs, allow_missing=True)[0]


def output_simjob_filename(config, **kwargs):
    """Returns the full path to the output file for a `simid`, `tier` and job index."""
    tier = kwargs.get("tier")

    if tier is None:
        msg = "the 'tier' argument is mandatory"
        raise RuntimeError(msg)

    fname = simjob_rel_basename() + f"-tier_{tier}" + FILETYPES["output"][tier]
    expr = str(Path(config["paths"][f"tier_{tier}"]) / fname)
    return expand(expr, **kwargs, allow_missing=True)[0]


def output_simjob_regex(config, **kwargs):
    tier = kwargs.get("tier")

    if tier is None:
        msg = "the 'tier' argument is mandatory"
        raise RuntimeError(msg)

    fname = "*-tier_{tier}" + FILETYPES["output"][tier]
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


def smk_ver_filename_for_stp(config, wildcards):
    """Returns the vertices file needed for the 'stp' tier job, if needed. Used
    as lambda function in the `build_tier_stp` Snakemake rule."""
    tdir = template_macro_dir(config, tier="stp")

    with (tdir / "simconfig.yaml").open() as f:
        sconfig = yaml.safe_load(f)[wildcards.simid]

    if "vertices" in sconfig:
        return output_simjob_filename(config, tier="ver", simid=sconfig["vertices"])
    else:
        return []


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


def pdf_config_path(config):
    return template_macro_dir(config, tier="pdf") / "build-pdf-config.yaml"


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
