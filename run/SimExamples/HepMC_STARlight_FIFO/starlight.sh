#! /usr/bin/env bash
#
# The script must only write HepMC to stdout

STARLIGHT_ROOT=$(starlight-config)
cp $STARLIGHT_ROOT/config/slight.in .
rm slight.out
# run STARlight (it needs slight.in)
starlight &> slight.out &
PROC1=$!
tail -f slight.out | awk -f ${STARLIGHT_ROOT}/HepMC/pdgMass.awk -f ${STARLIGHT_ROOT}/HepMC/starlight2hepmc.awk &
PROC2=$!

wait ${PROC1}
kill ${PROC2}
