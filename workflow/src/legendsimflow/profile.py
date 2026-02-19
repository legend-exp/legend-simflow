# Copyright (C) 2026 Luigi Pertoldi <gipert@pm.me>
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

import contextlib
import logging
import time
from collections import OrderedDict, defaultdict
from collections.abc import Callable

import psutil

log = logging.getLogger(__name__)


def _f(x: int | float) -> str:
    return f"{x:.5g}"


def _pct(x: int | float) -> str:
    return f"{x:.1f}"


def make_profiler() -> tuple[Callable, Callable]:
    proc = psutil.Process()
    stats = OrderedDict()
    stats = defaultdict(
        lambda: {
            "wall_s": 0.0,
            # per-call rss deltas (rss_after - rss_before);
            # we keep both the maximum delta and the average delta
            "max_delta_rss_mb": float("-inf"),
            "avg_delta_rss_mb": 0.0,
            "n_calls": 0,
        }
    )

    @contextlib.contextmanager
    def profile_block(name: str) -> None:
        t0 = time.perf_counter()
        rss0 = proc.memory_info().rss

        try:
            yield
        finally:
            t1 = time.perf_counter()
            rss1 = proc.memory_info().rss

            delta_mb = (rss1 - rss0) / 1024**2

            s = stats[name]
            s["wall_s"] += t1 - t0

            # update max delta
            s["max_delta_rss_mb"] = max(s["max_delta_rss_mb"], delta_mb)

            # online mean update for the average delta
            s["n_calls"] += 1
            s["avg_delta_rss_mb"] += (delta_mb - s["avg_delta_rss_mb"]) / s["n_calls"]

            msg = (
                f"block {name}: delta_rss_mb={_f(delta_mb)} wall_time_s={_f(t1 - t0)} "
            )
            log.debug(msg)

    def print_stats() -> None:
        msg = "==== profiling report ===="
        log.info(msg)

        blocks = list(stats)
        sum_wall_s = sum(stats[b]["wall_s"] for b in blocks)
        sum_max_delta_rss_mb = sum(stats[b]["max_delta_rss_mb"] for b in blocks)
        sum_avg_delta_rss_mb = sum(stats[b]["avg_delta_rss_mb"] for b in blocks)

        for block in blocks:
            s = stats[block]

            wall_s = s["wall_s"]
            max_delta_rss_mb = s["max_delta_rss_mb"]
            avg_delta_rss_mb = s["avg_delta_rss_mb"]

            wall_s_frac = 100.0 * wall_s / sum_wall_s if sum_wall_s else 0.0
            max_delta_rss_mb_frac = (
                100.0 * max_delta_rss_mb / sum_max_delta_rss_mb
                if sum_max_delta_rss_mb
                else 0.0
            )
            avg_delta_rss_mb_frac = (
                100.0 * avg_delta_rss_mb / sum_avg_delta_rss_mb
                if sum_avg_delta_rss_mb
                else 0.0
            )

            msg = (
                f"block {block} ]]] "
                f"wall_time_s={_f(wall_s)} ({_pct(wall_s_frac)}%) "
                f"max_delta_rss_mb={_f(max_delta_rss_mb)} ({_pct(max_delta_rss_mb_frac)}%) "
                f"avg_delta_rss_mb={_f(avg_delta_rss_mb)} ({_pct(avg_delta_rss_mb_frac)}%)"
            )
            log.info(msg)
        msg = "=========================="
        log.info(msg)

    return profile_block, print_stats
