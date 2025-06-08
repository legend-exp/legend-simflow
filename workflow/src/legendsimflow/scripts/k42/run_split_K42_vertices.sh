#!/bin/bash

doit() {
    if [ -f "generated/tier/raw/l200a-lar-inside-ms-K42/l200a-lar-inside-ms-K42_$1-tier_raw.root" ] && \
       [ -f "generated/tier/raw/l200a-lar-outside-ms-K42/l200a-lar-outside-ms-K42_$1-tier_raw.root" ]; then
        return
    fi
    echo "INFO: spawning job $1"
    shifter --image=legendexp/legend-software:latest -- root -l -q -x \
        "./workflow/scripts/k42/split_K42_vertices.C(
        \"generated/tier/raw/l200a-lar-cylinder-K42/l200a-lar-cylinder-K42_$1-tier_raw.root\",
        \"generated/tier/raw/l200a-lar-inside-ms-K42/l200a-lar-inside-ms-K42_$1-tier_raw.root\",
        \"generated/tier/raw/l200a-lar-outside-ms-K42/l200a-lar-outside-ms-K42_$1-tier_raw.root\")"
}
export -f doit

nfiles=$(find generated/tier/raw/l200a-lar-cylinder-K42 -maxdepth 1 -mindepth 1 | wc -l)
seq -f "%04g" 0 $((nfiles-1)) | parallel doit {}
