#!/usr/bin/env bash

macro="$1"
logfile="$2"

MaGe --no-colors "$macro" &> "$logfile"

code=0

grep -i "command not found" "$logfile" && code=1
grep -B 10 -i "batch is interrupted" "$logfile" && code=1
grep -B 10 -P "^(Error|Fatal):\w+" "$logfile" && code=1
grep -C 100 -iP "segmentation (fault|violation)" "$logfile" && code=1
grep -m 1 "you will get the same vertex position from now on!" "$logfile" && code=1
grep -A 15 -- "-------- EEEE ------- G4Exception-START -------- EEEE -------" "$logfile" && code=1
grep -C 13 -i "track stuck or not moving" "$logfile" && code=1

if [ "$code" != "0" ]; then
   echo -e "\n"
   echo -e " ┌──────────────────────────────────────────────────────────┐"
   echo -e " │ MaGe returned a non-zero code, see above or the log file │"
   echo -e " └──────────────────────────────────────────────────────────┘"
fi

exit $code
