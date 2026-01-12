# ruff: noqa: I002

# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>,
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

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
import legendmeta
import matplotlib.pyplot as plt

from legendsimflow import hpge_pars, nersc, utils
from legendsimflow import metadata as mutils
from legendsimflow.plot import decorate

args = snakemake  # nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

config = args.config

if "l200data" not in args.config.paths:
    msg = "Cannot extract current pars without setting the path to the l200data in the simconfig file."
    raise KeyError(msg)

l200data = args.config.paths.l200data

runid = args.wildcards.runid
hpge = args.wildcards.hpge_detector
metadata = args.config.metadata
pars_file = args.output.pars_file
plot_file = args.output.plot_file
log_file = args.log[0]

# setup logging
logger = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)
hit_tier_name = utils.get_hit_tier_name(l200data)

msg = f"... determined hit tier name is {hit_tier_name}"
logger.info(msg)
logger.info("... looking up the fit inputs")

raw_file, wf_idx, dsp_cfg_file = hpge_pars.lookup_currmod_fit_inputs(
    l200data,
    metadata,
    runid,
    hpge,
    hit_tier_name,
)

lh5_group = mutils._get_lh5_table(
    metadata,
    raw_file,
    hpge,
    "raw",
    runid,
)

logger.info("... fetching the current pulse")
t, A = hpge_pars.get_current_pulse(raw_file, lh5_group, wf_idx, str(dsp_cfg_file))

logger.info("... fitting the current pulse to extract the model")
popt, x, y = hpge_pars.fit_currmod(t, A)

# now plot
logger.info("... plotting the fit result")
fig, _ = hpge_pars.plot_currmod_fit_result(t, A, x, y)
decorate(fig)
plt.savefig(plot_file)

logger.info("... adding the mean aoe")
mean_aoe = hpge_pars.estimate_mean_aoe(popt)
popt_dict = utils._curve_fit_popt_to_dict(popt)

# logger.info("... estimating effect of noise")
# a_resolution = hpge_pars.

logger.info("... saving outputs")
dbetto.utils.write_dict(
    {"current_pulse_pars": popt_dict, "mean_aoe": mean_aoe}, pars_file
)
