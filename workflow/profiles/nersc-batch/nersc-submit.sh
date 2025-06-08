#!/bin/bash

function sbatch_submit() {
    local job="$1"
    local logdir="$2"
    shift 2

    if [[ -n "$DRY_RUN" ]]; then
        echo "INFO: WOULD HAVE SUBMITTED $job with args:"
        echo "INFO:     $*"
    return
    else
        echo "INFO: SUBMITTING $job with args:"
        echo "INFO:     $*"
    fi

    sbatch \
        --nodes 1 \
        --ntasks-per-node=1 \
        --account m2676 \
        --constraint cpu \
        --time 12:00:00 \
        --qos regular \
        --licenses scratch,cfs \
        --job-name "$job" \
        --output "$logdir/$job.log" \
        --error "$logdir/$job.log" \
        --image legendexp/legend-software:latest \
        --wrap "
            srun snakemake \
                --shadow-prefix $PSCRATCH \
                $*
        "
}

function sbatch_submit_simid() {
    local simid="$1"
    local version="$2"
    local logdir="$3"

    local job="$version::$simid"

    echo "INFO: inspecting $job"
    is_job_in_queue "$job" && return

    snakemake --config simlist="$simid" --dry-run 2> /dev/null | grep -q 'Nothing to be done' && return

    sbatch_submit "$job" "$logdir" --config simlist="$simid"
}

function is_job_in_queue() {
    if squeue --me --format '%200j' | grep -q "$1"; then
        echo "INFO: job already queued, skipping"
        return 0
    else return 1
    fi
}

export -f sbatch_submit
export -f sbatch_submit_simid
export -f is_job_in_queue

version="$(basename "$PWD")"
version="$(basename "$(dirname "$PWD")")-$version"

logdir=".slurm/$(date +'%Y%m%dT%H%M%SZ')"
if [[ -z "$DRY_RUN" ]]; then
    mkdir -p "$logdir"
    ln -fns "$(basename "$logdir")" .slurm/latest
fi

if [[ "$1" == "parallel" ]]; then
    simids=$(python -c '
import json

with open("inputs/simprod/config/tier/stp/l200a/simconfig.json") as f:
    simids = json.load(f).keys()

for s in simids:
    print(f"pdf.{s}", end=" ")
    ')

    # shellcheck disable=SC2086
    parallel sbatch_submit_simid ::: $simids ::: "$version" ::: "$logdir"

else
    job="$version::all"

    echo "INFO: inspecting $job"
    is_job_in_queue "$job" && exit 1

    if [[ -z "$NO_CHECK" ]]; then
        snakemake --dry-run "$@" 2> /dev/null | grep -q 'Nothing to be done' && exit 1
    fi

    sbatch_submit "$job" "$logdir" "$*"
fi

echo "INFO: logs in $(realpath "$logdir")"
