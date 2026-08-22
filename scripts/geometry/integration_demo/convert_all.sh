#!/bin/bash
# Stage 1 of the integration demo: convert IRIS and Bagger, each twice.
#
#   <model>_exact  : the cascade  CSG -> exact O2BVHSurfaceSolid -> tessellated fallback
#   <model>_tess   : pure tessellation, same mesh precision as the cascade's fallback
#
# Usage: convert_all.sh <output-root>
set -u
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}
export PYTHONPATH=${PYTHONPATH:-}
OUT=${1:?usage: convert_all.sh <output-root>}
GEO=$(cd "$(dirname "$0")/.." && pwd)
source "$GEO/integration_demo/env_converter.sh"

# Mesh precision. --mesh-prec sets linear AND angular deflection to the same value and
# Stream_P_RepresentationBench.md measured it behaving as an *angular* knob. The default 0.1
# is unaffordable on IRIS (one 2 m sphere reached 22.9 GB resident and was killed);
# 0.5 is the finest value measured to fit in memory on the full model.
IRIS_PREC=${IRIS_PREC:-0.5}
BAGGER_PREC=${BAGGER_PREC:-0.1}

NIST="$GEO/g4_nist_database/G4_NIST_DB.json"
run() {  # run <tag> <args...>
  local tag=$1; shift
  mkdir -p "$OUT/conv/$tag" "$OUT/logs"
  echo "=== $tag ==="
  /usr/bin/time -v $PYOCC "$GEO/O2_CADtoTGeo.py" "$@" \
      --output-folder "$OUT/conv/$tag" -o geom.C --g4-nist-json "$NIST" \
      > "$OUT/logs/conv_$tag.log" 2>&1
  echo "  exit=$? -> $OUT/logs/conv_$tag.log"
}

run bagger_exact "$GEO/STEP_examples/Bagger.step" --csg auto --exact-surfaces auto \
    --materials-csv "$GEO/integration_demo/Bagger_MATERIALS.csv"
run bagger_tess  "$GEO/STEP_examples/Bagger.step" --mesh --mesh-prec "$BAGGER_PREC" \
    --materials-csv "$GEO/integration_demo/Bagger_MATERIALS.csv"
run iris_exact   "$GEO/ALICE_3_example/CAD_noETA.stp" --csg auto --exact-surfaces auto \
    --mesh --mesh-prec "$IRIS_PREC" --materials-csv "$GEO/IRIS/IRIS_MATERIALS.csv"
run iris_tess    "$GEO/ALICE_3_example/CAD_noETA.stp" --mesh --mesh-prec "$IRIS_PREC" \
    --materials-csv "$GEO/IRIS/IRIS_MATERIALS.csv"


# The exact-surface macros need one post-processing step before o2-sim can JIT them; see
# patch_exact_macro.py for the precise reason (a converter/loader namespace interaction).
python3 "$GEO/integration_demo/patch_exact_macro.py" "$OUT"/conv/*/geom.C
echo "done"
