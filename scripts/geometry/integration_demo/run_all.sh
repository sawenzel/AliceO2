#!/bin/bash
# The measurement matrix of the integration demo. Two geometry representations of the same
# two CAD models, everything else held fixed: same seed, same generator, single-threaded,
# MCStepLogger attached (per-step ROOT tree included).
#
#   geantino_*  : 5 x 50 geantinos, eta [-1,1] -- pure geometry, no interactions
#   electron_*  : 5 x 20 electrons at 1 GeV
#   pion_*      : 5 x 20 pi+ at 1 GeV
#   matfan_*    : one event of NRAYS geantinos on Fibonacci-sphere directions, the
#                 per-ray material-budget equivalence test
#
# Usage: run_all.sh <conv-root>
set -u
GEO=$(cd "$(dirname "$0")" && pwd)
OUT=${1:?usage: run_all.sh <conv-root>}
EV=${EV:-5}
N=${N:-50}
NCHG=${NCHG:-20}
NRAYS=${NRAYS:-512}

for rep in exact tess; do
  EVENTS=$EV NGUN=$N    PDG=0   STEPLOG=1 "$GEO/run_sim.sh" "$OUT" $rep geantino_$rep
  EVENTS=$EV NGUN=$NCHG PDG=11  STEPLOG=1 "$GEO/run_sim.sh" "$OUT" $rep electron_$rep
  EVENTS=$EV NGUN=$NCHG PDG=211 STEPLOG=1 "$GEO/run_sim.sh" "$OUT" $rep pion_$rep
  EVENTS=1 GEN=extgen STEPLOG=1 \
    CONFIGKEY="GeneratorExternal.fileName=$GEO/fibonacci_geantinos.macro;GeneratorExternal.funcName=fibonacci_geantinos($NRAYS,1.0,0.9)" \
    "$GEO/run_sim.sh" "$OUT" $rep matfan_$rep
done
echo "run_all done"
