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
from __future__ import annotations

import argparse
import logging
import math
from datetime import timedelta

from dbetto import utils
from dbetto.time import str_to_datetime
from legendmeta import LegendMetadata, LegendSlowControlDB


def round_step_5(x):
    return int(5 * math.floor(x / 5 + 0.5))


parser = argparse.ArgumentParser()

parser.add_argument("runsel", help="run selection string, i.e. l200-p14-r001-phy")
parser.add_argument("-o", "--output", help="output file path", required=True)

args = parser.parse_args()


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)

log.info("initializing...")

experiment, period, run, datatype = args.runsel.split("-")

# need to handle db names, damn ivano...
db_name = "scdbL140" if int(period[1:]) < 13 else "scdb"

lmeta = LegendMetadata()
scdb = LegendSlowControlDB()
scdb.connect(db_name=db_name)

timestamp = lmeta.datasets.runinfo[period][run][datatype].start_key
msg = f"start timestamp of {args.runsel} is {timestamp}"
log.info(msg)

# add 30 minutes to avoid any slow control time lag
timestamp = str_to_datetime(timestamp) + timedelta(minutes=30)

chmap = lmeta.channelmap(timestamp).group("system").geds.map("name")

log.info("querying the LEGEND Slow Control...")

voltages = {}
for name in sorted(chmap):
    meta = chmap[name]
    status = scdb.status(meta, on=timestamp)

    msg = f"voltage set/mon for channel {name} is {status.vset}/{status.vmon} V"
    logging.info(msg)

    if abs(status.vset - status.vmon) > 5:
        logging.warning("set and monitored voltage differ by more than 5 V!")

    # use monitored voltage: I noticed that sometimes vset can be different
    # from vmon for a long time period (and vmon is always correct)
    voltages[name] = {"operational_voltage_in_V": round_step_5(status.vmon)}

utils.write_dict(voltages, args.output)
