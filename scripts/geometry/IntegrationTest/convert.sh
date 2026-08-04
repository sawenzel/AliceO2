#!/bin/bash
# Stage 1+2 of the CAD->TGeo Geant integration test (Handoff_IntegrationTest.md):
# four conversions -- {ALICE3 IRIS, oTOF} x {exact cascade, pure tessellated} -- with
# materials assigned (ALICE3 from the IRIS BOM, oTOF one material for the whole module).
#
# Usage:   convert.sh <workdir> [which]
#   <workdir>  output root; each conversion goes to its own subfolder
#   [which]    optional: one of alice3_exact alice3_mesh otof_exact otof_mesh (default: all)
#
# Environment expected (see Handoff_IntegrationTest.md section 11): alienv O2 dev environment
# loaded, LD_LIBRARY_PATH with $B/stage first then OCCT/Python prepended, PYTHONPATH with
# pythonOCC prepended, O2_BUILD_DIR set, and the alibuild python3.10 on PATH as $PYOCC.
#
# Mesh precision: 0.5 everywhere meshing happens.
#   - NOT the default 0.1: meshing ALICE3 at 0.1 was measured at 22.9 GB resident on one
#     2 m sphere and killed (NEXT.md); 0.5 and 1.0 were measured affordable (679-841 MB).
#     0.5 is the finer of the two affordable values.
#   - Same value for both representations, so the exact run's tessellated-fallback parts are
#     identical to the pure-mesh run's parts and the comparison isolates the 20 exact solids.
#   - The precision *sweep* is a roadmap item, deliberately not this task.
#
# The exact configurations also pass --mesh: the converter cascade is CSG -> exact surfaces
# -> tessellated, and without --mesh the fallback for a part that emits neither CSG nor
# sidecar is a bounding-box placeholder (12 triangles), not the part. ALICE3 has 35 such
# solids; a physics comparison against that would be meaningless.

set -euo pipefail

WORK=${1:?usage: convert.sh <workdir> [which]}
WHICH=${2:-all}
GEODIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYOCC=${PYOCC:-python3}
MESHPREC=0.5

ALICE3_STEP="$GEODIR/ALICE_3_example/CAD_noETA.stp"
OTOF_STEP="$GEODIR/STEP_examples/oTOF System V3-R92cm.step"
ALICE3_MAT=(--materials-csv "$GEODIR/IRIS/IRIS_MATERIALS.csv" --g4-nist-json "$GEODIR/g4_nist_database/G4_NIST_DB.json")
OTOF_MAT=(--materials-csv "$GEODIR/IntegrationTest/otof_materials.csv" --g4-nist-json "$GEODIR/g4_nist_database/G4_NIST_DB.json")

EXACT=(--csg auto --exact-surfaces auto --mesh --mesh-prec $MESHPREC --dump-brep)
MESH=(--mesh --mesh-prec $MESHPREC)

run() {
  local name=$1; shift
  local step=$1; shift
  if [[ "$WHICH" != "all" && "$WHICH" != "$name" ]]; then return 0; fi
  local out="$WORK/$name"
  mkdir -p "$out"
  echo "=== $name -> $out"
  ( cd "$GEODIR/.." && /usr/bin/time -v $PYOCC "$GEODIR/O2_CADtoTGeo.py" "$step" \
      --output-folder "$out" -o geom.C \
      --surface-report "$out/surface_report.json" \
      "$@" ) > "$out/convert.log" 2>&1
  echo "=== $name done, exit $?"
}

run alice3_exact "$ALICE3_STEP" "${EXACT[@]}" "${ALICE3_MAT[@]}"
run alice3_mesh  "$ALICE3_STEP" "${MESH[@]}"  "${ALICE3_MAT[@]}"
run otof_exact   "$OTOF_STEP"   "${EXACT[@]}" "${OTOF_MAT[@]}"
run otof_mesh    "$OTOF_STEP"   "${MESH[@]}"  "${OTOF_MAT[@]}"
echo "all requested conversions finished"
