from __future__ import annotations

import argparse
import random
import signal
import subprocess
from pathlib import Path

import yaml
from dbetto import AttrsDict
from legenddataflowscripts.workflow.utils import subst_vars
from legendmeta import LegendMetadata

from . import aggregate


def _partition(xs, n):
    k, r = divmod(len(xs), n)
    out, i = [], 0
    for j in range(n):
        s = k + (j < r)
        out.append(xs[i : i + s])
        i += s
    return out


def snakemake_nersc_cli():
    parser = argparse.ArgumentParser(
        description="Execute the Simflow on multiple nodes in parallel."
    )
    parser.add_argument(
        "-N", "--nodes", type=int, required=True, help="number of nodes"
    )
    parser.add_argument(
        "--without-srun",
        action="store_true",
        help="do not prefix the snakemake call with 'srun ...'",
    )
    args, extra = parser.parse_known_args()

    if args.nodes < 2:
        msg = "must parallelize over at least 2 nodes"
        raise ValueError(msg)

    cfg_path = Path("./simflow-config.yaml")
    if not cfg_path.is_file():
        msg = "this program must be executed in the directory where simflow-config.yaml resides"
        raise RuntimeError(msg)

    with cfg_path.open("r") as f:
        config = yaml.safe_load(f)

    subst_vars(
        config,
        var_values={"_": Path().resolve()},
        use_env=True,
        ignore_missing=False,
    )
    config = AttrsDict(config)

    # NOTE: this will attempt a clone of legend-metadata, if the directory does not exist
    metadata = LegendMetadata(config.paths.metadata, lazy=True)

    if "legend_metadata_version" in config:
        metadata.checkout(config.legend_metadata_version)

    config["metadata"] = metadata

    simlist = config.get("simlist", None)
    make_tiers = config.make_tiers
    if simlist is None:
        # auto determine tier from config
        tiers = ("pdf", "cvt", "evt", "hit", "opt", "stp")
        tier = next(t for t in tiers if t in make_tiers)

        simlist = [
            f"{tier}.{simid}" for simid in aggregate.gen_list_of_all_simids(config)
        ]

    # trick: there won't be anything to do for some simids (targets already
    # done), this could result in a very inefficient partitioning. as a
    # mitigation, we randomly shuffle the simlist first
    random.shuffle(simlist)

    procs = []
    for simlist_chunk in _partition(simlist, args.nodes):
        smk_cmd = [
            "snakemake",
            "--workflow-profile",
            "workflow/profiles/nersc",
            "--config",
            "simlist=" + ",".join(simlist_chunk),
            *extra,
        ]
        if not args.without_srun:
            smk_cmd = [
                "srun",
                "--disable-status",  # otherwise SIGINT has no effect
                "--nodes",
                "1",
                "--ntasks",
                "1",
                "--cpus-per-task",
                "256",
                *smk_cmd,
            ]

        print("INFO: spawning process:", " ".join(smk_cmd))  # noqa: T201
        procs.append(subprocess.Popen(smk_cmd))

    # propagate signals to the snakemake instances.
    def new_signal_handler(sig: int, _):
        for p in procs:
            p.send_signal(sig)

    signals = [
        signal.SIGHUP,
        signal.SIGINT,
        signal.SIGQUIT,
        signal.SIGTERM,
        signal.SIGTSTP,  # SIGSTOP cannot be caught, and will do nothing...
        signal.SIGCONT,
        signal.SIGUSR1,
        signal.SIGUSR2,
        signal.SIGWINCH,
    ]

    for sig in signals:
        signal.signal(sig, new_signal_handler)

    for p in procs:
        rc = p.wait()
        if rc != 0:
            msg = f"process failed: {p.args}"
            raise RuntimeError(msg)

    print("INFO: all snakemake processes successfully returned")  # noqa: T201
