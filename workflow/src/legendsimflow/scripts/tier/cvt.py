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

import shutil

import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
from lgdo import lh5

from legendsimflow import nersc

args = nersc.dvs_ro_snakemake(snakemake)  # noqa: F821

evt_files = args.input
cvt_file = args.output[0]
log_file = args.log[0]
metadata = args.config.metadata

BUFFER_LEN = "500*MB"


# setup logging
log = ldfs.utils.build_log(metadata.simprod.config.logging, log_file)

if len(evt_files) == 1:
    shutil.copy(evt_files[0], cvt_file)
else:
    for table in lh5.ls(evt_files[0]):
        for chunk in lh5.LH5Iterator(evt_files, table, buffer_len=BUFFER_LEN):
            lh5.write(chunk, table, cvt_file, wo_mode="append")
