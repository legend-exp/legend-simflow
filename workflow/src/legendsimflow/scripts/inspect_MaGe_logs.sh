#!/usr/bin/env bash

# shellcheck disable=SC2154
logs="${snakemake_params[logdir]}"

function it() {
    code=0

    grep -A 15 -i -- "-------- WWWW ------- G4Exception-START -------- WWWW -------" "$1" && code=1
    grep -P "^Warning:\w+" "$1" && code=1

    if [ "$code" != "0" ]; then
        echo
        echo -e "\033[1m^^^^^ $1 ^^^^^\033[0m"
        echo
    fi
}

find "$logs" -regextype "egrep" -regex ".*(ver|raw)/.+/.+.log" \
    | while read -r file; do it "$file"; done
