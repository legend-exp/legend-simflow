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

import sys
from unittest.mock import patch

import pytest
import yaml

from legendsimflow.cli import (
    _partition,
    snakemake_nersc_batch_cli,
    snakemake_nersc_cli,
)

# ---------------------------------------------------------------------------
# _partition
# ---------------------------------------------------------------------------


def test_partition_even():
    result = _partition(list(range(6)), 3)
    assert result == [[0, 1], [2, 3], [4, 5]]


def test_partition_uneven():
    result = _partition(list(range(7)), 3)
    # 7 // 3 = 2 remainder 1 → first chunk gets 3 items
    assert [len(c) for c in result] == [3, 2, 2]
    assert sum(len(c) for c in result) == 7
    # All original items present exactly once
    assert sorted(x for chunk in result for x in chunk) == list(range(7))


def test_partition_single_chunk():
    items = [1, 2, 3]
    result = _partition(items, 1)
    assert result == [items]


def test_partition_empty():
    assert _partition([], 3) == [[], [], []]


# ---------------------------------------------------------------------------
# snakemake_nersc_cli
# ---------------------------------------------------------------------------


def _make_nersc_config(tmp_path, simlist=("stp.simA", "stp.simB")):
    """Write a minimal simflow-config.yaml suitable for snakemake_nersc_cli tests."""
    cfg = {
        "paths": {"metadata": str(tmp_path)},
        "simlist": list(simlist),
    }
    (tmp_path / "simflow-config.yaml").write_text(yaml.dump(cfg))


def test_nersc_cli_raises_for_single_node(monkeypatch):
    """ValueError is raised before the config-file check, so no config is needed."""
    monkeypatch.setattr(sys, "argv", ["snakemake-nersc", "-N", "1"])
    with pytest.raises(ValueError, match="at least 2 nodes"):
        snakemake_nersc_cli()


def test_nersc_cli_raises_missing_config(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["snakemake-nersc", "-N", "2"])
    with pytest.raises(RuntimeError, match=r"simflow-config\.yaml"):
        snakemake_nersc_cli()


def test_nersc_cli_no_submit(tmp_path, monkeypatch, capsys):
    """--no-submit prints the snakemake commands that would be spawned."""
    _make_nersc_config(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["snakemake-nersc", "-N", "2", "--no-submit"])
    with patch("legendsimflow.cli.LegendMetadata"):
        snakemake_nersc_cli()
    out = capsys.readouterr().out
    assert "snakemake" in out
    # default: srun prefix is included
    assert "srun" in out


def test_nersc_cli_without_srun(tmp_path, monkeypatch, capsys):
    """--without-srun omits the srun prefix from each command."""
    _make_nersc_config(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        ["snakemake-nersc", "-N", "2", "--no-submit", "--without-srun"],
    )
    with patch("legendsimflow.cli.LegendMetadata"):
        snakemake_nersc_cli()
    out = capsys.readouterr().out
    assert "snakemake" in out
    assert "srun" not in out


def test_nersc_cli_extra_args_forwarded(tmp_path, monkeypatch, capsys):
    """Extra arguments are forwarded verbatim to every snakemake invocation."""
    _make_nersc_config(tmp_path, simlist=["stp.simA"])
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        ["snakemake-nersc", "-N", "2", "--no-submit", "--dryrun", "--cores", "4"],
    )
    with patch("legendsimflow.cli.LegendMetadata"):
        snakemake_nersc_cli()
    out = capsys.readouterr().out
    assert "--dryrun" in out
    assert "--cores" in out


def test_nersc_cli_simlist_partitioned(tmp_path, monkeypatch, capsys):
    """The simlist is split across the requested number of nodes."""
    _make_nersc_config(tmp_path, simlist=["stp.s1", "stp.s2", "stp.s3", "stp.s4"])
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys, "argv", ["snakemake-nersc", "-N", "2", "--no-submit", "--without-srun"]
    )
    with patch("legendsimflow.cli.LegendMetadata"):
        snakemake_nersc_cli()
    out = capsys.readouterr().out
    # Two "would spawn" lines - one per node
    assert out.count("would spawn") == 2


# ---------------------------------------------------------------------------
# snakemake_nersc_batch_cli
# ---------------------------------------------------------------------------


def test_batch_cli_raises_missing_config(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["snakemake-nersc-batch", "-t", "02:00:00"])
    with pytest.raises(RuntimeError, match=r"simflow-config\.yaml"):
        snakemake_nersc_batch_cli()


def test_batch_cli_no_submit_single_node(tmp_path, monkeypatch, capsys):
    """Single node: uses plain snakemake (not snakemake-nersc)."""
    (tmp_path / "simflow-config.yaml").write_text("{}")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys, "argv", ["snakemake-nersc-batch", "-t", "02:00:00", "--no-submit"]
    )
    snakemake_nersc_batch_cli()
    out = capsys.readouterr().out
    assert "sbatch" in out
    assert "02:00:00" in out
    # single node → uses plain snakemake, not snakemake-nersc
    assert "snakemake-nersc" not in out


def test_batch_cli_no_submit_multi_node(tmp_path, monkeypatch, capsys):
    """Multiple nodes: uses snakemake-nersc."""
    (tmp_path / "simflow-config.yaml").write_text("{}")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        ["snakemake-nersc-batch", "-t", "02:00:00", "-N", "4", "--no-submit"],
    )
    snakemake_nersc_batch_cli()
    out = capsys.readouterr().out
    assert "sbatch" in out
    assert "snakemake-nersc" in out


def test_batch_cli_no_submit_with_job_name(tmp_path, monkeypatch, capsys):
    """-J/--job-name is included in the sbatch command."""
    (tmp_path / "simflow-config.yaml").write_text("{}")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        ["snakemake-nersc-batch", "-t", "02:00:00", "-J", "my-job", "--no-submit"],
    )
    snakemake_nersc_batch_cli()
    out = capsys.readouterr().out
    assert "my-job" in out


def test_batch_cli_no_submit_with_mail(tmp_path, monkeypatch, capsys):
    """--mail-user adds --mail-type and the address to the sbatch command."""
    (tmp_path / "simflow-config.yaml").write_text("{}")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "snakemake-nersc-batch",
            "-t",
            "02:00:00",
            "--mail-user",
            "user@example.com",
            "--no-submit",
        ],
    )
    snakemake_nersc_batch_cli()
    out = capsys.readouterr().out
    assert "user@example.com" in out
    assert "--mail-type" in out


def test_batch_cli_no_submit_cpus_per_task(tmp_path, monkeypatch, capsys):
    """-c/--cpus-per-task is passed through to sbatch."""
    (tmp_path / "simflow-config.yaml").write_text("{}")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "snakemake-nersc-batch",
            "-t",
            "02:00:00",
            "-c",
            "128",
            "--no-submit",
        ],
    )
    snakemake_nersc_batch_cli()
    out = capsys.readouterr().out
    assert "128" in out


def test_batch_cli_no_submit_extra_args_forwarded(tmp_path, monkeypatch, capsys):
    """Extra arguments beyond the known flags are forwarded to snakemake."""
    (tmp_path / "simflow-config.yaml").write_text("{}")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "snakemake-nersc-batch",
            "-t",
            "02:00:00",
            "--no-submit",
            "--dryrun",
        ],
    )
    snakemake_nersc_batch_cli()
    out = capsys.readouterr().out
    assert "--dryrun" in out
