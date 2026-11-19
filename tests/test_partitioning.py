from __future__ import annotations

import pytest

from legendsimflow.partitioning import partition_simstat, primaries_per_run


def test_partition_simstat_expected_dict():
    n_events = {"job_000": 10, "job_001": 10}
    n_events_part = {
        "l200-p03-r001-phy": 5,
        "l200-p03-r002-phy": 10,
        "l200-p03-r003-phy": 5,
    }
    runlist = list(n_events_part.keys())

    expected = {
        "job_000": {
            "l200-p03-r001-phy": [0, 4],
            "l200-p03-r002-phy": [5, 9],
        },
        "job_001": {
            "l200-p03-r002-phy": [0, 4],
            "l200-p03-r003-phy": [5, 9],
        },
    }

    assert partition_simstat(n_events, n_events_part, runlist) == expected


def test_partition_simstat_unsorted_runlist():
    n_events = {"job_000": 10, "job_001": 10}
    n_events_part = {
        "l200-p03-r001-phy": 5,
        "l200-p03-r002-phy": 10,
        "l200-p03-r003-phy": 5,
    }
    sorted_runlist = sorted(n_events_part.keys())
    reversed_runlist = list(reversed(sorted_runlist))

    assert partition_simstat(
        n_events, n_events_part, sorted_runlist
    ) == partition_simstat(n_events, n_events_part, reversed_runlist)


def test_partition_simstat_rejects_duplicate_runids():
    n_events = {"job_000": 10, "job_001": 10}
    n_events_part = {
        "l200-p03-r001-phy": 5,
        "l200-p03-r002-phy": 10,
        "l200-p03-r003-phy": 5,
    }
    # a duplicate would be visited twice and silently overwrite the first chunk
    runlist = [*n_events_part.keys(), "l200-p03-r002-phy"]

    with pytest.raises(ValueError, match="duplicate runids"):
        partition_simstat(n_events, n_events_part, runlist)


def test_primaries_per_run():
    # 10 / 3 each, the leftover primary goes to the first run
    livetimes = {"r1": 1.0, "r2": 1.0, "r3": 1.0}
    assert primaries_per_run(livetimes, 10) == {"r1": 4, "r2": 3, "r3": 3}
    assert primaries_per_run({"r1": 300.0, "r2": 100.0}, 1000) == {"r1": 750, "r2": 250}
    assert primaries_per_run({"r1": 1.0, "r2": 2.0}, 0) == {"r1": 0, "r2": 0}
    with pytest.raises(ValueError, match="zero"):
        primaries_per_run({"r1": 0.0}, 10)
