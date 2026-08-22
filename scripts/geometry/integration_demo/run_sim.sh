#!/bin/bash
# One o2-sim run of the integration demo.
#
# Usage: run_sim.sh <conv-root> <exact|tess> <run-tag> [extra o2-sim args...]
#
# Environment knobs:
#   EVENTS=3  SEED=42  GEN=boxgen  PDG=0 (geantino)  NGUN=20  STEPLOG=0|1  NOGEANT=0|1
#   CONFIGKEY="a=1;b=2"   extra --configKeyValues, appended to the box-gun ones
#
# Everything is deterministic: the seed is fixed and the run is single-threaded
# (o2-sim-serial), so two runs differing only in the geometry representation are comparable.
set -u
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}
GEO=$(cd "$(dirname "$0")/.." && pwd)
CONV=${1:?usage: run_sim.sh <conv-root> <exact|tess> <run-tag> [args...]}
REP=${2:?}
TAG=${3:?}
shift 3

EVENTS=${EVENTS:-3}
SEED=${SEED:-42}
GEN=${GEN:-boxgen}
PDG=${PDG:-0}
NGUN=${NGUN:-20}
PMIN=${PMIN:-1.0}
PMAX=${PMAX:-1.0}
STEPLOG=${STEPLOG:-0}
CONFIGKEY=${CONFIGKEY:-}
NOGEANT=${NOGEANT:-0}

RUNDIR=$CONV/runs/$TAG
mkdir -p "$RUNDIR"
python3 "$GEO/integration_demo/make_configs.py" "$CONV" "$REP" "$RUNDIR" || exit 1

source "$GEO/integration_demo/env_o2.sh" >/dev/null 2>&1
cd "$RUNDIR" || exit 1

ARGS=(-n "$EVENTS" -g "$GEN" --seed "$SEED"
      --detectorList "EXTCAD:$RUNDIR/detectorlist.json"
      --extGeomFile "$RUNDIR/externalGeometry.json"
      --configKeyValues "BoxGun.number=$NGUN;BoxGun.pdg=$PDG;BoxGun.prange[0]=$PMIN;BoxGun.prange[1]=$PMAX${CONFIGKEY:+;$CONFIGKEY}"
      -o o2sim)
[ "$NOGEANT" = "1" ] && ARGS+=(--noGeant)

if [ "$STEPLOG" = "1" ]; then
  export LD_PRELOAD=$SW/MCStepLogger/latest/lib/libMCStepLoggerInterceptSteps.so
  export MCSTEPLOG_OUTFILE=$RUNDIR/MCStepLoggerOutput.root
  # the per-step ROOT tree (StepLoggerTree) is only written when MCSTEPLOG_TTREE is set;
  # without it MCStepLogger only prints its per-volume summary to the log.
  export MCSTEPLOG_TTREE=1
fi

/usr/bin/time -v "$B/stage/bin/o2-sim-serial" "${ARGS[@]}" "$@" > "$RUNDIR/sim.log" 2>&1
rc=$?
unset LD_PRELOAD
echo "run $TAG ($REP) exit=$rc -> $RUNDIR/sim.log"
exit $rc
