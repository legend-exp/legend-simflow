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

# NOTE: this script is not ready to be integrated in Snakemake, since it only
# works on the LNGS cluster! At the end of this file there's an example of how
# it should look like

# ruff: noqa: F821, T201


import json
import math
from pathlib import Path

from legendmeta import JsonDB, LegendMetadata

prod_dir = "/data2/public/prodenv/prod-blind/ref-v1.0.0"
lmeta = LegendMetadata(f"{prod_dir}/inputs")

runlist = (
    [f"l200-p03-r00{r}-phy" for r in range(6)]
    + [f"l200-p04-r00{r}-phy" for r in range(4)]
    + [f"l200-p06-r00{r}-phy" for r in range(6)]
    + [f"l200-p07-r00{r}-phy" for r in range(1, 8)]
    + [f"l200-p08-r00{r}-phy" for r in range(5)]
    + [f"l200-p08-r00{r}-phy" for r in range(6, 10)]
)

# use FCCD values reviewed by Elisabetta
with Path("fccd-reviewed.json").open() as f:
    fccd = json.load(f)["fccd-mm"]

# get generated parameters (for energy resolution) from the data production
par_pht_meta = JsonDB(f"{prod_dir}/generated/par/pht", lazy=True)

for run in runlist:
    print(">>>", run)

    # get parameters and hardware configuration for run
    p = run.split("-")
    tstamp = lmeta.dataprod.runinfo[p[1]][p[2]][p[3]].start_key
    chmap = lmeta.channelmap(tstamp).group("system").geds
    hit_pars = par_pht_meta.on(tstamp)

    # now loop over all HPGe channels
    evt_cfg = {}
    for _, data in chmap.items():
        # compute the MaGe sensitive volume identifier
        mage_id = int(1010000 + 100 * data.location.string + data.location.position)

        # get energy resolution curves
        eres_pars = [None, None, None]
        channel = f"ch{data.daq.rawid}"
        if channel in hit_pars:
            pars = None
            # first try getting curves specific to the run
            # i.e. parameters of sqrt(a + b*E)
            try:
                ene_cfg = hit_pars[channel].results.ecal.cuspEmax_ctc_cal
                pars = [*list(ene_cfg.eres_linear.parameters.values()), 0]
            except AttributeError:
                # if not found, use curves specific to partition
                try:
                    pars = hit_pars[
                        channel
                    ].results.partition_ecal.cuspEmax_ctc_cal.eres_linear.parameters.values()
                except AttributeError:
                    if data.analysis.usability not in ("off", "ac"):
                        print(
                            f"WARNING: no eres params found for {data.analysis.usability} {channel}/{data.name} in {run}"
                        )

            # convert to format expected by mage-post-proc
            if pars is not None and all(p >= 0 for p in pars):
                eres_pars = [round(math.sqrt(x) / 2.355, 6) for x in pars]
        elif data.analysis.usability not in ("off", "ac"):
            print(
                f"ERROR: {data.analysis.usability} {channel}/{data.name} not in JSON file"
            )

        # finally add JSON block about channel
        evt_cfg[mage_id] = {
            "name": data.name,
            "nplus-fccd-mm": fccd[data.name],
            "energy": dict(zip(["sig0", "sig1", "sig2"], eres_pars)),
            "usability": data.analysis.usability,
        }

    with Path(f"{run}-build_evt.json").open("w") as f:
        json.dump(evt_cfg, f, indent=2)


# import os
# from pathlib import Path

# # TODO: generate the files from metadata inplace here
# # this script is run for each data taking run
# # snakemake.input[0] points to the legend-metadata clone in legend-simflow
# # snakemake.output[0] is the output filename

# src = (
#     Path(snakemake.input[0])
#     / "simprod"
#     / "config"
#     / "tier"
#     / "evt"
#     / snakemake.config["experiment"]
#     / f"{snakemake.wildcards.runid}-build_evt.json"
# )
# dest = snakemake.output[0]

# os.system(f"cp {src} {dest}")
