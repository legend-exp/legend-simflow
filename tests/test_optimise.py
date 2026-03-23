from __future__ import annotations

import math

import pytest

from legendsimflow import optimise


def test_round_down_nsf_large():
    assert optimise.round_down_nsf(1291860, 2) == 1200000
    assert optimise.round_down_nsf(717.7, 2) == 710
    assert optimise.round_down_nsf(1800000, 2) == 1800000
    assert optimise.round_down_nsf(1801000, 2) == 1800000
    assert optimise.round_down_nsf(100, 2) == 100
    assert optimise.round_down_nsf(999, 2) == 990


def test_round_down_nsf_small():
    # values with <= n sig figs: just truncate
    assert optimise.round_down_nsf(3.6, 2) == 3
    assert optimise.round_down_nsf(9.9, 2) == 9
    assert optimise.round_down_nsf(31, 2) == 31
    assert optimise.round_down_nsf(99, 2) == 99


def test_round_down_nsf_zero_and_negative():
    assert optimise.round_down_nsf(0, 2) == 0
    assert optimise.round_down_nsf(-5, 2) == 0


def test_round_up_nsf_large():
    assert optimise.round_up_nsf(314, 2) == 320
    assert optimise.round_up_nsf(3714, 2) == 3800
    assert optimise.round_up_nsf(3700, 2) == 3700
    assert optimise.round_up_nsf(100, 2) == 100
    assert optimise.round_up_nsf(101, 2) == 110


def test_round_up_nsf_small():
    # values with <= n sig figs: just take ceiling
    assert optimise.round_up_nsf(3, 2) == 3
    assert optimise.round_up_nsf(3.7, 2) == 4
    assert optimise.round_up_nsf(31, 2) == 31
    assert optimise.round_up_nsf(99, 2) == 99


def test_round_up_nsf_zero_and_negative():
    assert optimise.round_up_nsf(0, 2) == 0
    assert optimise.round_up_nsf(-5, 2) == 0


def test_round_down_rounds_down():
    # round_down should never increase the value
    for x in [1.9, 50.0, 99.9, 100.1, 999.0, 9999.9, 100000.0]:
        result = optimise.round_down_nsf(x, 2)
        assert result <= x, f"round_down_nsf({x}) = {result} > {x}"


def test_round_up_rounds_up():
    # round_up should never decrease the value
    for x in [1.0, 50.0, 99.0, 100.1, 999.0, 9999.9, 100001.0]:
        result = optimise.round_up_nsf(x, 2)
        assert result >= x, f"round_up_nsf({x}) = {result} < {x}"


def test_round_nsf_significant_figures():
    # results should have at most n significant figures
    # i.e., all digits below the (n-1)-th significant position are zero
    def has_at_most_n_sig_figs(x: int, n: int) -> bool:
        if x == 0:
            return True
        k = math.floor(math.log10(x))
        trailing_zeros_needed = max(0, k - n + 1)
        return x % (10**trailing_zeros_needed) == 0

    for x in [314, 3714, 12345, 99999, 100001]:
        result = optimise.round_up_nsf(x, 2)
        assert has_at_most_n_sig_figs(result, 2), (
            f"round_up_nsf({x}, 2) = {result} has more than 2 sig figs"
        )

    for x in [314, 3714, 12345, 99999, 100001]:
        result = optimise.round_down_nsf(x, 2)
        assert has_at_most_n_sig_figs(result, 2), (
            f"round_down_nsf({x}, 2) = {result} has more than 2 sig figs"
        )


def test_compute_primaries_per_job_basic():
    # speed=1000 ev/s, max_runtime=30min -> raw=1800000 -> round_down -> 1800000
    result = optimise.compute_primaries_per_job(1000, 30, n_sig_figs=2)
    assert result == 1800000

    # speed=500 ev/s, max_runtime=30min -> raw=900000 -> round_down -> 900000
    result = optimise.compute_primaries_per_job(500, 30, n_sig_figs=2)
    assert result == 900000

    # result must be <= raw (no job exceeds time limit)
    speed, max_rt = 717.7, 30
    raw = speed * max_rt * 60
    result = optimise.compute_primaries_per_job(speed, max_rt, n_sig_figs=2)
    assert result <= raw


def test_compute_primaries_per_job_constraint():
    # Job must not exceed max_runtime_mins
    for speed in [10.0, 100.0, 1000.0, 5000.0]:
        for max_rt in [10, 30, 60]:
            ppj = optimise.compute_primaries_per_job(speed, max_rt, n_sig_figs=2)
            actual_runtime_mins = ppj / speed / 60
            assert actual_runtime_mins <= max_rt, (
                f"speed={speed}, max_rt={max_rt}: runtime {actual_runtime_mins:.1f} > {max_rt}"
            )


def test_compute_number_of_jobs_basic():
    # ceil(1_000_000 / 900_000) = 2, round_up(2, 2) = 2
    result = optimise.compute_number_of_jobs(900_000, 1_000_000, n_sig_figs=2)
    assert result == 2

    # ceil(1_000_000 / 100_000) = 10, round_up(10, 2) = 10
    result = optimise.compute_number_of_jobs(100_000, 1_000_000, n_sig_figs=2)
    assert result == 10

    # ceil(1_000_001 / 1_000_000) = 2, round_up(2, 2) = 2
    result = optimise.compute_number_of_jobs(1_000_000, 1_000_001, n_sig_figs=2)
    assert result == 2


def test_compute_number_of_jobs_constraint():
    # Total primaries must be >= target
    for ppj in [1000, 10_000, 100_000, 1_000_000]:
        for target in [1_000_000, 10_000_000, 100_000_000]:
            njobs = optimise.compute_number_of_jobs(ppj, target, n_sig_figs=2)
            assert ppj * njobs >= target, (
                f"ppj={ppj}, target={target}: {ppj}*{njobs}={ppj * njobs} < {target}"
            )


def test_compute_sim_settings_basic():
    ppj, njobs = optimise.compute_sim_settings(
        speed=500, max_runtime_mins=30, target_total_primaries=1_000_000
    )
    assert isinstance(ppj, int)
    assert isinstance(njobs, int)
    assert ppj > 0
    assert njobs > 0


def test_compute_sim_settings_constraints():
    speed = 500.0  # ev/s
    max_rt = 30  # min
    target = 10_000_000

    ppj, njobs = optimise.compute_sim_settings(speed, max_rt, target)

    # Constraint 1: no job exceeds max runtime
    actual_runtime_mins = ppj / speed / 60
    assert actual_runtime_mins <= max_rt

    # Constraint 2: total primaries >= target
    assert ppj * njobs >= target


def test_compute_sim_settings_minimises_jobs():
    speed = 1000.0
    max_rt = 30
    target = 5_000_000

    ppj, njobs = optimise.compute_sim_settings(speed, max_rt, target)

    # Verify we cannot reduce njobs by 1 and still meet target
    # (unless njobs is already 1)
    if njobs > 1:
        assert ppj * (njobs - 1) < target


def test_compute_sim_settings_docstring_example():
    ppj, njobs = optimise.compute_sim_settings(
        speed=500, max_runtime_mins=30, target_total_primaries=1_000_000
    )
    assert ppj == 900000
    assert njobs == 2


@pytest.mark.parametrize(
    ("speed", "max_rt", "target"),
    [
        (0.1, 30, 100),
        (10.0, 60, 100_000),
        (10_000.0, 10, 1_000_000_000),
    ],
)
def test_compute_sim_settings_various(speed, max_rt, target):
    ppj, njobs = optimise.compute_sim_settings(speed, max_rt, target)
    assert ppj >= 1
    assert njobs >= 1
    assert ppj / speed / 60 <= max_rt
    assert ppj * njobs >= target
