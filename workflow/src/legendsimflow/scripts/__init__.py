# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>
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
import shlex


def log_script_invocation(
    log: logging.Logger,
    pixi_task: str,
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> None:
    """Log the equivalent CLI command to reproduce the current script invocation.

    Reconstructs the command from the parsed argument namespace so it works
    correctly both when called standalone and when run via Snakemake (where
    ``sys.argv`` reflects the Snakemake runner, not the script arguments).

    Parameters
    ----------
    log
        Logger to write the message to.
    pixi_task
        Name of the pixi task (e.g. ``"tier-cvt"``).
    parser
        The script's argument parser, used to map destinations back to flags.
    args
        The parsed argument namespace.
    """
    parts = [f"pixi run {pixi_task}"]
    for action in parser._actions:
        if not action.option_strings:
            continue
        val = getattr(args, action.dest, None)
        if val is None:
            continue
        flag = action.option_strings[0]
        if isinstance(val, list):
            parts += [flag] + [shlex.quote(str(v)) for v in val]
        elif isinstance(val, bool):
            if val:
                parts.append(flag)
        else:
            parts += [flag, shlex.quote(str(val))]
    log.info("CLI invocation: %s", " ".join(parts))
