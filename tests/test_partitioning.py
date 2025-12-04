from legendsimflow.partitioning import partition_simstat


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
