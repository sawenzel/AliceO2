#!/bin/bash
# Stage 4 of the CAD->TGeo Geant integration test: the simulation runs, one representation
# at a time, everything seeded and few-event so a full run is minutes.
#
# Usage: run_sim.sh <workdir> <mode> <repr>
#   <workdir>  the workdir used by convert.sh / make_sim_configs.sh
#   <mode>     geomcheck | pions | steplog | matbudget
#   <repr>     exact | mesh
#
# Modes:
#   geomcheck  o2-sim-serial -n 0: build + close the full assembled geometry, transport
#              nothing (the Stage-3 "check before simulating").
#   pions      5 events x 10 pi (BoxGun defaults: pdg 211, eta [-1,1], p 0.1-5 GeV),
#              TGeant4 with geomRoot (TGeoNavigator, the representation actually under test).
#   steplog    the same run under MCStepLogger (LD_PRELOAD interception) for per-volume
#              step tallies.
#   matbudget  o2-sim-evalmat, one event containing a geantino scan over an eta-phi grid
#              (matm.neta x matm.nphi = 60 x 360 directions, NOT an axis raster -- a
#              parallel/axis beam is direction-poor and has bitten this project twice).
#              Writes o2sim_matbudget.root with X0 / lambda / g-cm2 maps.
#
# Environment: alienv O2 dev env + LD_LIBRARY_PATH with $B/stage FIRST (the branch's
# freshly built binaries and libraries; `which o2-sim` resolves to the stale install
# otherwise) + VMCWORKDIR pointing at this source tree + ROOT_INCLUDE_PATH carrying
# Detectors/Base/include (the exact geom.C includes branch headers the installed O2 does
# not have). Do NOT carry the converter's pythonOCC/PYTHONPATH prepends into these runs;
# they poison the embedded Python and o2-sim segfaults at startup.

set -euo pipefail
WORK=${1:?usage: run_sim.sh <workdir> <mode> <repr>}
MODE=${2:?mode: geomcheck|pions|steplog|matbudget}
REPR=${3:?repr: exact|mesh}

: "${B:?source the sim environment first (B, stage-first LD_LIBRARY_PATH, VMCWORKDIR)}"
SEED=424242
NEV=${NEV:-5}   # events for pions/steplog; override via environment
COMMON=(--detectorList "INTEG:$WORK/sim/detlist.json" -m A3IRIS OTOF
        --extGeomFile "$WORK/sim/extgeom_$REPR.json" --seed $SEED)

OUT="$WORK/sim/${MODE}_${REPR}"
mkdir -p "$OUT"
cd "$OUT"

case "$MODE" in
  geomcheck)
    /usr/bin/time -v "$B/stage/bin/o2-sim-serial" -n 0 -e TGeant4 -g boxgen \
      "${COMMON[@]}" -o "$MODE" > run.log 2>&1
    ;;
  pions)
    /usr/bin/time -v "$B/stage/bin/o2-sim-serial" -n $NEV -e TGeant4 -g boxgen \
      "${COMMON[@]}" -o "$MODE" > run.log 2>&1
    ;;
  steplog)
    MCSL="$HOME/alisw/sw/ubuntu2404_aarch64/MCStepLogger/latest/lib/libMCStepLoggerInterceptSteps.so"
    [ -f "$MCSL" ] || { echo "MCStepLogger not found: $MCSL"; exit 3; }
    LD_PRELOAD="$MCSL" /usr/bin/time -v "$B/stage/bin/o2-sim-serial" -n $NEV -e TGeant4 -g boxgen \
      "${COMMON[@]}" -o "$MODE" > run.log 2>&1
    ;;
  matbudget)
    /usr/bin/time -v "$B/stage/bin/o2-sim-evalmat" -n 1 -e TGeant4 -g boxgen \
      "${COMMON[@]}" \
      --configKeyValues "matm.neta=60;matm.etamin=-1.5;matm.etamax=1.5;matm.nphi=360;matm.rmax=200;matm.zmax=500" \
      -o "$MODE" > run.log 2>&1
    ;;
  *) echo "unknown mode $MODE"; exit 2 ;;
esac
echo "$MODE/$REPR finished, exit $? -- output in $OUT"
