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

# ruff: noqa: F821, T201

import json
import re
from pathlib import Path

import snakemake as smk
from utils import aggregate, patterns, utils

with Path(snakemake.input.cfgfile).open() as f:
    config = json.load(f)[snakemake.params.simid]

n_prim = config.pop("number_of_primaries")
n_macros = None
outver_list = None
is_benchmark = False

# determine whether external vertices are required
if "vertices" in config:
    # get list of ver output files
    outver_list = aggregate.gen_list_of_simid_outputs(
        snakemake.config, "ver", config.pop("vertices")
    )

    # if number of jobs is not specified, use number of vertices files
    if "number_of_jobs" not in config:
        n_macros = len(outver_list)

# if n_macros could not be determined, there MUST be an explicit reference to
# it in the config
if not n_macros:
    n_macros = config.pop("number_of_jobs")

if "benchmark" in snakemake.config:
    is_benchmark = snakemake.config["benchmark"].get("enabled", False)
    if is_benchmark:
        n_macros = 1
        n_prim = snakemake.config["benchmark"]["n_primaries"][snakemake.params.tier]

# prepare global substitution rules
substitutions = {"NUMBER_OF_PRIMARIES": int(n_prim / n_macros)}

# now that we processed the "special" configuration fields, we interpret all
# the others as variables to be expanded in the template macro
for k, v in config.items():
    if isinstance(v, str):
        substitutions.update({k.upper(): v})
    elif isinstance(v, list) and all(isinstance(s, str) for s in v):
        substitutions.update({k.upper(): "\n".join(v)})
    else:
        print(f"WARNING: key '{k}' does not define a valid substitution rule")

# first substitute global variables
with Path(snakemake.input.template).open() as f:
    text = utils.subst_vars(f.read().strip(), substitutions, ignore_missing=True)

    # if benchmark run, write all events to disk, not only those that deposit
    # energy in the detectors. TODO: hardcoding application-specific commands
    # is not nice here
    if is_benchmark:
        text = re.sub(
            r"\n */MG/io/MCRun/setWriteEventsThatDepositEnergyInGe +true *\n",
            "\n",
            text,
        )

# then substitute macro-specific variables
for i in range(n_macros):
    # determine output file name for this macro
    outname = patterns.output_simjob_filename(
        snakemake.config,
        tier=snakemake.params.tier,
        simid=snakemake.params.simid,
        jobid=f"{i:>04d}",
    )

    substitutions.update({"OUTPUT_FILE": outname})
    # might also need to set a vertices file name (see above)
    if outver_list:
        substitutions.update({"VERTICES_FILE": outver_list[i]})

    # determine the macro file name for write out
    inname = Path(
        patterns.input_simjob_filename(
            snakemake.config,
            tier=snakemake.params.tier,
            simid=snakemake.params.simid,
            jobid=f"{i:>04d}",
        )
    )

    text_out = utils.subst_vars(text, substitutions).strip()

    # check if the file exists and open it
    # NOTE: unfortunately, this doesn't produce the desired effect with
    # Snakemake, as it always removes output files before regenerating them:
    # https://stackoverflow.com/questions/49178033/snakemake-avoid-removing-output-files-before-executing-the-shell-command
    if inname.is_file():
        with inname.open() as f:
            # if file content is already the expected one, do not overwrite it
            if f.read().strip() == text_out:
                continue

    # otherwise, prepare for writing and write
    smk.utils.makedirs(str(inname.parent))
    with inname.open("w") as f:
        f.write(text_out)
