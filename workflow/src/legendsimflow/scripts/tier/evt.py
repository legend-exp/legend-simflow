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

from pathlib import Path

import legenddataflowscripts as ldfs
import legenddataflowscripts.utils

args = snakemake  # noqa: F821

opt_file = args.input.opt_file
hit_file = args.input.hit_file
evt_file = args.output[0]
log_file = args.log[0]
metadata = args.config.metadata
buffer_len = args.params.buffer_len


# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)

Path(evt_file).touch()
