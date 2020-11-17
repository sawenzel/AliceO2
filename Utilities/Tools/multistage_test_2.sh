# a simple test for the development of GRID checkpointing
# and mixing 1-core with multi-core hardware within the same job workflow

. ${O2_ROOT}/share/scripts/jobutils.sh

JOBUTILS_SKIPDONE=ON

echo "MY ALIEN WORKDIR BELONGING TO THIS JOB = ${ALIEN_JOB_OUTPUTDIR}"

# THIS IS THE SIMULATION TASK
taskwrapper simstage.log "o2-sim-serial -m ITS TPC -g boxgen -n 100"

# THIS IS THE DIGITIZATION TASK
for i in `seq 0 20`; do
  taskwrapper digistage_${i}.log "o2-sim-digitizer-workflow -b --onlyDet ITS"
done

