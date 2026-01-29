from __future__ import annotations

import argparse
import random
import shlex
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
    """Implementation of the ``snakemake-nersc`` CLI."""
    parser = argparse.ArgumentParser(
        description="""Execute the Simflow on multiple nodes in parallel.
        Extra arguments will be forwarded to Snakemake."""
    )
    parser.add_argument(
        "--no-submit",
        action="store_true",
        help="do not run anything, just show what would be run.",
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

    simlist = list(config.get("simlist", None))  # make a copy for shuffling later
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
            "workflow/profiles/nersc-compute",
            "--config",
            "simlist=" + ",".join(simlist_chunk),
            *extra,
        ]
        if not args.without_srun:
            smk_cmd = [
                "srun",
                "--nodes",
                "1",
                "--ntasks",
                "1",
                "--cpus-per-task",
                "256",
                *smk_cmd,
            ]

        if not args.no_submit:
            print("INFO: spawning process:", shlex.join(smk_cmd))  # noqa: T201
            procs.append(subprocess.Popen(smk_cmd))
        else:
            print("INFO: would spawn:", " ".join(smk_cmd))  # noqa: T201

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


def snakemake_nersc_batch_cli():
    """Implementation of the ``snakemake-nersc-batch`` CLI."""
    parser = argparse.ArgumentParser(
        description="""Execute the Simflow as a batch Slurm job.
        See sbatch docs for help about the CLI arguments.
        Extra arguments will be forwarded to Snakemake."""
    )
    parser.add_argument(
        "--no-submit",
        action="store_true",
        help="do not run sbatch, just show what would be run.",
    )

    parser.add_argument("-t", "--time", required=True)
    parser.add_argument("-N", "--nodes", default="1")
    parser.add_argument("-c", "--cpus-per-task", default="256")
    parser.add_argument("-J", "--job-name")
    parser.add_argument("--mail-user")

    args, smk_args = parser.parse_known_args()

    cfg_path = Path("./simflow-config.yaml")
    if not cfg_path.is_file():
        msg = "this program must be executed in the directory where simflow-config.yaml resides"
        raise RuntimeError(msg)

    cmd = [
        "sbatch",
        "--qos",
        "regular",
        "--constraint",
        "cpu",
        "--ntasks",
        args.nodes,
        "--ntasks-per-node",
        "1",
        "--account",
        "m2676",
        "--licenses",
        "cfs,scratch",
        "--output",
        ".slurm/jobid-%j.log",
        "--error",
        ".slurm/jobid-%j.log",
        "--time",
        args.time,
        "--nodes",
        args.nodes,
        "--cpus-per-task",
        args.cpus_per_task,
    ]
    if args.job_name is not None:
        cmd.extend(["--job-name", args.job_name])
    if args.mail_user is not None:
        cmd.extend(["--mail-type", "TIME_LIMIT,FAIL", "--mail-user", args.mail_user])

    if int(args.nodes) > 1:
        snakemake = "pixi run snakemake-nersc --nodes $SLURM_NNODES"
    else:
        snakemake = (
            "pixi run snakemake --workflow-profile workflow/profiles/nersc-compute"
        )

    cmd.extend(
        [
            "--wrap",
            f"{snakemake} --keep-going " + " ".join(smk_args),
        ]
    )

    if args.no_submit:
        print("INFO: would run:", shlex.join(cmd))  # noqa: T201

    else:
        print("INFO: running:", shlex.join(cmd))  # noqa: T201
        rc = subprocess.Popen(cmd).wait()
        if rc != 0:
            msg = "submission failed"
            raise RuntimeError(msg)
