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
from pathlib import Path

from dbetto import TextDB, utils
from dbetto.time import datetime_to_str, str_to_datetime
from legendmeta import LegendMetadata, LegendSlowControlDB


def round_step_5(x):
    return int(5 * math.floor(x / 5 + 0.5))


def dict_diff(d1: dict, d2: dict) -> dict:
    """return entries of d1 that differ from d2"""
    return {k: v for k, v in d1.items() if k not in d2 or d2[k] != v}


parser = argparse.ArgumentParser()

parser.add_argument("runsel", help="run selection string, i.e. l200-p14-r001-phy")
parser.add_argument("-o", "--output", help="output file path")
parser.add_argument("--opv-db", help="path to existing opv database (folder)")

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
    log.info(msg)

    if abs(status.vset - status.vmon) > 5:
        log.warning("set and monitored voltage differ by more than 5 V!")

    # use monitored voltage: I noticed that sometimes vset can be different
    # from vmon for a long time period (and vmon is always correct)
    voltages[name] = {"operational_voltage_in_V": round_step_5(status.vmon)}

# now we check what we need to write
if args.opv_db is not None:
    # load the voltages in the database for the current timestamp
    prev_voltages = TextDB(args.opv_db).on(timestamp)
    # compute the diff
    v_diff = dict_diff(voltages, prev_voltages)

    if v_diff:
        msg = "voltages have changed compared to existing values, updating the database"
        log.info(msg)

        # write down the diff only
        opv_filename = f"{experiment}-{period}-{run}-T%-all-opvs.yaml"
        utils.write_dict(v_diff, Path(args.opv_db) / opv_filename)

        # and add validity entry
        validity = utils.load_dict(Path(args.opv_db) / "validity.yaml")
        validity.append(
            {
                "valid_from": datetime_to_str(timestamp),
                "mode": "append",
                "apply": [opv_filename],
            }
        )
        # sort the validity file
        validity = sorted(validity, key=lambda x: x["valid_from"])
        utils.write_dict(validity, Path(args.opv_db) / "validity.yaml")
    else:
        msg = "voltages have not changed compared to existing values! continuing"
        log.info(msg)

else:
    utils.write_dict(voltages, args.output)
