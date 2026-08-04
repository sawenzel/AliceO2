#!/bin/bash
# Stage 5 helper: side-by-side per-volume step/secondary tallies from two MCStepLogger runs.
# Usage: summarize_steplog.sh <exact_run.log> <mesh_run.log>
# Parses the [STEPLOGGER] per-event 'VolName <name> COUNT <n> SECONDARIES <m>' lines
# (volume names may contain spaces) and aggregates over events.
set -euo pipefail
A=${1:?usage: summarize_steplog.sh <exact_run.log> <mesh_run.log>}
B=${2:?}

parse() {
  grep -E "VolName .* COUNT" "$1" \
    | sed 's/.*VolName //; s/ COUNT / | /; s/ SECONDARIES / | /; s/ P\[.*//' \
    | awk -F' \\| ' '{c[$1]+=$2; s[$1]+=$3} END {for (v in c) printf "%s|%d|%d\n", v, c[v], s[v]}'
}

join -t'|' -a1 -a2 -e 0 -o 0,1.2,1.3,2.2,2.3 \
  <(parse "$A" | sort -t'|' -k1,1) <(parse "$B" | sort -t'|' -k1,1) \
  | awk -F'|' 'BEGIN {printf "%-30s %10s %8s %10s %8s\n", "volume", "stepsA", "secA", "stepsB", "secB"}
               {printf "%-30s %10d %8d %10d %8d\n", $1, $2, $3, $4, $5; ta+=$2; tb+=$4; sa+=$3; sb+=$5}
               END {printf "%-30s %10d %8d %10d %8d\n", "TOTAL", ta, sa, tb, sb}'
echo
for f in "$A" "$B"; do
  echo "$f:"
  grep -oE "did [0-9]+ steps" "$f" | awk '{a+=$2} END {print "  steps (logger total):", a}'
  grep -oE "transported [0-9]+ different tracks" "$f" | awk '{a+=$2} END {print "  tracks:", a}'
done
