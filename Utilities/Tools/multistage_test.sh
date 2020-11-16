# a simple test for the development of GRID checkpointing
# and mixing 1-core with multi-core hardware within the same job workflow

. ${O2_ROOT}/share/scripts/jobutils.sh

JOBUTILS_SKIPDONE=ON

echo "MY ALIEN WORKDIR BELONGING TO THIS JOB = ${ALIEN_JOB_OUTPUTDIR}"

# THIS IS THE SIMULATION TASK
taskwrapper simstage.log "o2-sim-serial -m ITS TPC -g boxgen"

# this is code that should finally be "injected by the driver" based on an annotation (or function call)
# pragma ALIEN-CHECKPOINT SITE
if [ "${ALIEN_JOB_OUTPUTDIR}" ]; then
  taskwrapper checkpoint1.log "touch checkpoint1.log_done ; tar -czf checkpoint.tar.gz * ; alien.py cp -f checkpoint.tar.gz alien://${ALIEN_JOB_OUTPUTDIR}/; exit 0"
fi

# THIS IS THE DIGITIZATIOn TASK
taskwrapper digistage.log "o2-sim-digitizer-workflow -b --onlyDet ITS"

# we should do some cleanup here ---> do this in driver?
if [ -f "digistage.log" ]; then
  alien.py rm alien://${ALIEN_JOB_OUTPUTDIR}/checkpoint.tar.gz
fi
