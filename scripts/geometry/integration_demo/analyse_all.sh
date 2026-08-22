#!/bin/bash
# Reduce every run's MCStepLogger tree to a text tally, and every geantino run to a per-ray
# material budget. Writes <conv-root>/analysis/{steps_<tag>.txt,matbudget_<tag>.txt}.
#
# Usage: analyse_all.sh <conv-root>
set -u
GEO=$(cd "$(dirname "$0")" && pwd)
OUT=${1:?usage: analyse_all.sh <conv-root>}
mkdir -p "$OUT/analysis"
source "$GEO/env_o2.sh" >/dev/null 2>&1
# MCStepLogger is not on the O2 environment's paths; only the analysis needs it, and it is
# deliberately kept out of env_o2.sh so the sim runs see exactly what they saw when measured.
export LD_LIBRARY_PATH=$SW/MCStepLogger/latest/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$SW/MCStepLogger/latest/include:$ROOT_INCLUDE_PATH
for d in "$OUT"/runs/*/; do
  tag=$(basename "$d")
  sf="$d/MCStepLoggerOutput.root"
  [ -f "$sf" ] || continue
  root -l -b -q "$GEO/analyse_steps.C(\"$sf\")" > "$OUT/analysis/steps_$tag.txt" 2>&1
  case "$tag" in
    geantino_*|matfan_*)
      root -l -b -q "$GEO/matbudget.C(\"$d/o2sim_geometry.root\",\"$sf\",\"$OUT/analysis/matbudget_$tag.txt\")" \
        > "$OUT/analysis/matbudget_$tag.log" 2>&1
      ;;
  esac
  echo "analysed $tag"
done
