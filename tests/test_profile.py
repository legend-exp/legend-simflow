from __future__ import annotations

from legendsimflow.profile import make_profiler


def test_make_profiler_returns_three_callables():
    profile_block, print_stats, print_stats_since_last = make_profiler()

    assert callable(profile_block)
    assert callable(print_stats)
    assert callable(print_stats_since_last)


def test_profiler_context_manager_and_printers_are_callable():
    profile_block, print_stats, print_stats_since_last = make_profiler()

    with profile_block("block_a"):
        pass

    with profile_block("block_b"):
        pass

    print_stats()
    print_stats_since_last()
    print_stats_since_last()
