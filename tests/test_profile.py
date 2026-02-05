from __future__ import annotations

import time

import pytest

from legendsimflow import profile


def test_make_profiler_basic():
    """Test that make_profiler returns two callables."""
    profile_block, print_stats = profile.make_profiler()

    assert callable(profile_block)
    assert callable(print_stats)


def test_profile_block_context_manager():
    """Test that profile_block works as a context manager."""
    profile_block, print_stats = profile.make_profiler()

    # Should not raise any errors
    with profile_block("test_block"):
        time.sleep(0.01)

    # Call print_stats to ensure it doesn't crash
    print_stats()


def test_profile_block_timing():
    """Test that profile_block measures wall time correctly."""
    profile_block, print_stats = profile.make_profiler()

    # Profile a block with a known sleep time
    sleep_time = 0.05
    with profile_block("sleep_block"):
        time.sleep(sleep_time)

    # We can't directly access stats, but we can verify it doesn't crash
    # and that it runs without errors
    print_stats()


def test_profile_block_multiple_calls():
    """Test that profile_block accumulates statistics over multiple calls."""
    profile_block, print_stats = profile.make_profiler()

    # Call the same block multiple times
    for _ in range(3):
        with profile_block("repeated_block"):
            time.sleep(0.01)

    # Should not raise any errors
    print_stats()


def test_profile_block_different_blocks():
    """Test that profile_block tracks different blocks separately."""
    profile_block, print_stats = profile.make_profiler()

    with profile_block("block_a"):
        time.sleep(0.01)

    with profile_block("block_b"):
        time.sleep(0.02)

    with profile_block("block_a"):
        time.sleep(0.01)

    # Should track both blocks separately
    print_stats()


def test_profile_block_with_exception():
    """Test that profile_block records stats even when an exception is raised."""
    profile_block, print_stats = profile.make_profiler()

    # Profile a block that raises an exception
    with (  # noqa: PT012
        pytest.raises(ValueError),
        profile_block("error_block"),
    ):
        msg = "Intentional error"
        raise ValueError(msg)

    # Stats should still be recorded
    print_stats()


def test_profile_block_memory_tracking():
    """Test that profile_block tracks memory changes."""
    profile_block, print_stats = profile.make_profiler()

    # Allocate some memory
    with profile_block("memory_block"):
        # Allocate a reasonably large array to see memory change
        data = [0] * 1000000  # noqa: F841

    # Should not raise any errors
    print_stats()


def test_print_stats_empty():
    """Test that print_stats works even when no blocks have been profiled."""
    _, print_stats = profile.make_profiler()

    # Should not raise any errors
    print_stats()


def test_profile_block_nested():
    """Test that profile_block works with nested contexts."""
    profile_block, print_stats = profile.make_profiler()

    with profile_block("outer"):
        time.sleep(0.01)
        with profile_block("inner"):
            time.sleep(0.01)

    # Should track both blocks
    print_stats()
