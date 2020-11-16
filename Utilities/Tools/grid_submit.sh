#!/bin/bash

# A helper script, making it easy to submit existing
# scripts as an ALIEN GRID job (with the following notation):
#
# grid-submit my_script.sh jobname
#
# The script then handles all interaction with the GRID automatically. The user
# does not need to create JDLs files nor upload them to the GRID manually/herself.
#
# The script can also simulate execution of the job locally. To this end, it suffices
# to say
#
# ALIEN_PROC_ID=1 grid-submit my_script.sh
#
# Currently handles only a very basic JDL configuration. Further improvements would be:
#
# -) allow JDL customization via command line arguments or JDL tags inside the script
#

# set -o pipefail

function per() { printf "\033[31m$1\033[m\n" >&2; }
function pok() { printf "\033[32m$1\033[m\n" >&2; }
function banner() { echo ; echo ==================== $1 ==================== ; }

function Usage() { echo "$0 --script scriptname | -c WORKDIR_RELATIVE_TO_TOP [ --jobname JOBNAME ] [ --topworkdir WORKDIR (ON TOP OF HOME) ] "; }

# find out if this script is really executed on GRID
# in this case, we should find an environment variable JALIEN_TOKEN_CERT
ONGRID=0
[ "${JALIEN_TOKEN_CERT}" ] && ONGRID=1

# All is redirected to log.txt but kept on stdout as well
if [[ $ALIEN_PROC_ID ]]; then
  exec &> >(tee -a alien_log_${ALIEN_PROC_ID}.txt)
fi

# this tells us to continue an existing job --> in this case we don't create a new workdir
while [ $# -gt 0 ] ; do
    case $1 in
	-c) CONTINUE_WORKDIR=$2;  shift 2 ;;   # this should be the workdir of a job to continue (without HOME and ALIEN_TOPWORKDIR)
        --local) LOCAL_MODE="ON"; shift 1 ;;   # if we want emulate execution in the local workdir (no GRID interaction)
        --script) SCRIPT=$2; shift 2 ;;  # the job script to submit
        --jobname) JOBNAME=$2; shift 2 ;; # the job name associated to the job --> determined directory name on GRID
        --topworkdir) ALIEN_TOPWORKDIR=$2; shift 2 ;; # the top workdir relative to GRID home
        --ttl) JOBTTL=$2; shift 2 ;; # allows to specifiy ttl for job
        --partition) GRIDPARTITION=$2; shift 2 ;; # allows to specificy a GRID partition for the job
        --dry) DRYRUN="ON"; shift 1 ;; # do a try run and not actually interact with the GRID (just produce local jdl file)
        --o2tag) O2TAG=$2; shift 2 ;; # 
	-h) Usage ; exit ;;
        *) break ;;
    esac
done

echo "SCRIPT ${SCRIPT}"

# analyse options:
# we should either run with --script or with -c
[ "${SCRIPT}" ] && [ "$CONTINUE_WORKDIR" ] && echo "Script and continue mode not possible at same time" && exit 1
if [ "${ONGRID}" = 0 ]; then
  [[ ! ( "${SCRIPT}" || "$CONTINUE_WORKDIR" ) ]] && echo "Either script or continue mode required" && exit 1
fi

# General job configuration
MY_USER=${ALIEN_USER:-`whoami`}
if [[ ! $MY_USER ]]; then
  per "Problems retrieving current AliEn user. Did you run alien-token-init?"
  exit 1
fi
MY_HOMEDIR="/alice/cern.ch/user/${MY_USER:0:1}/${MY_USER}"
MY_JOBPREFIX="$MY_HOMEDIR/${ALIEN_TOPWORKDIR:-selfjobs}"
MY_JOBSCRIPT="$(cd "$(dirname "${SCRIPT}")" && pwd -P)/$(basename "${SCRIPT}")" # the job script with full path
MY_JOBNAME=${JOBNAME:-$(basename ${MY_JOBSCRIPT})}
MY_JOBNAMEDATE="${MY_JOBNAME}-$(date -u +%Y%m%d-%H%M%S)"
MY_JOBWORKDIR="${MY_JOBPREFIX}/${MY_JOBNAMEDATE}"  # ISO-8601 UTC
[ "${CONTINUE_WORKDIR}" ] && MY_JOBWORKDIR="${MY_JOBPREFIX}/${CONTINUE_WORKDIR}"
MY_BINDIR="$MY_JOBWORKDIR"

pok "Your job's working directory will be $MY_JOBWORKDIR"
pok "Set the job name by running $0 <scriptname> <jobname>"

#
# Generate local workdir
#
if [[ "${ONGRID}" == "0" ]]; then
  WORKDIR=${WORKDIR:-/tmp/alien_work/$(basename "$MY_JOBWORKDIR")}
  [ ! -d "${WORKDIR}" ] && mkdir -p ${WORKDIR}
  [ ! "${CONTINUE_WORKDIR}" ] && cp "${MY_JOBSCRIPT}" "${WORKDIR}/alien_jobscript.sh"
fi

# 
# Submitter code (we need to submit whenever a script is given as input and we are not in local mode)
#
[[ ( ! "${LOCAL_MODE}" ) && ( "${SCRIPT}" || "${CONTINUE_WORKDIR}" ) ]] && IS_ALIEN_JOB_SUBMITTER=ON

if [[ "${IS_ALIEN_JOB_SUBMITTER}" ]]; then
  #  --> test if alien is there?
  which alien.py 2> /dev/null
  # check exit code
  if [[ ! "$?" == "0"  ]]; then
    XJALIEN_LATEST=`find /cvmfs/alice.cern.ch/el7-x86_64/Modules/modulefiles/xjalienfs -type f -printf "%f\n" | tail -n1`
    banner "Loading xjalienfs package $XJALIEN_LATEST since not yet loaded"
    eval "$(/cvmfs/alice.cern.ch/bin/alienv printenv xjalienfs::"$XJALIEN_LATEST")"
  fi

  # Create temporary workdir to assemble files, and submit from there (or execute locally)
  cd "$(dirname "$0")"
  THIS_SCRIPT="$PWD/$(basename "$0")"

  cd "${WORKDIR}"

  # ---- Generate JDL ----------------
  # TODO: Make this configurable or read from a preamble section in the jobfile
  cat > "${MY_JOBNAMEDATE}.jdl" <<EOF
Executable = "${MY_BINDIR}/${MY_JOBNAMEDATE}.sh";
Arguments = "${CONTINUE_WORKDIR:+"-c ${CONTINUE_WORKDIR}"} --local ${O2TAG:+--o2tag ${O2TAG}}";
InputFile = "LF:${MY_JOBWORKDIR}/alien_jobscript.sh";
OutputDir = "${MY_JOBWORKDIR}";
Output = {
  "logs*.zip,alien*.txt,*digit*.root@disk=2"
};
Requirements = member(other.GridPartitions,"${GRIDPARTITION:-cc7}");
MemorySize = "60GB";
TTL=${JOBTTL:-600};
EOF
#

  pok "Local working directory is $PWD"
  if [ ! "${DRYRUN}" ]; then
    command_file="alien_commands.txt"

    pok "Preparing job \"$MY_JOBNAMEDATE\""
    (
      set -x
      # assemble all GRID interaction in a single script / transaction
      [ -f "${command_file}" ] && rm ${command_file}
      [ ! "${CONTINUE_WORKDIR}" ] && echo "rmdir ${MY_JOBWORKDIR}" >> ${command_file}    # remove existing job dir
      # echo "mkdir ${MY_BINDIR}" >> ${command_file}                      # create bindir
      echo "mkdir ${MY_JOBPREFIX}" >> ${command_file}                   # create job output prefix
      [ ! "${CONTINUE_WORKDIR}" ] && echo "mkdir ${MY_JOBWORKDIR}" >> ${command_file}
      echo "rm ${MY_BINDIR}/${MY_JOBNAMEDATE}.sh" >> ${command_file}    # remove current job script
      echo "cp ${PWD}/${MY_JOBNAMEDATE}.jdl alien://${MY_JOBWORKDIR}/${MY_JOBNAMEDATE}.jdl" >> ${command_file}  # copy the jdl
      echo "cp ${THIS_SCRIPT} alien://${MY_BINDIR}/${MY_JOBNAMEDATE}.sh" >> ${command_file}  # copy current job script to AliEn
      [ ! "${CONTINUE_WORKDIR}" ] && echo "cp ${MY_JOBSCRIPT} alien://${MY_JOBWORKDIR}/alien_jobscript.sh" >> ${command_file}

#      [ ! "${CONTINUE_WORKDIR}" ] && alien.py rmdir "$MY_JOBWORKDIR" || true    # remove existing job dir
#      alien.py mkdir "$MY_BINDIR" || true                                       # create bindir
#      alien.py mkdir "$MY_JOBPREFIX" || true                                    # create job output prefix
#      # alien.py mkdir jdl || true
#      [ ! "${CONTINUE_WORKDIR}" ] && alien.py mkdir "$MY_JOBWORKDIR" || true
#      alien.py rm "$MY_BINDIR/${MY_JOBNAMEDATE}.sh" || true                     # remove current job script
#      alien.py cp "${PWD}/${MY_JOBNAMEDATE}.jdl" alien://${MY_JOBWORKDIR}/${MY_JOBNAMEDATE}.jdl@ALICE::CERN::EOS || true  # copy the jdl
#      alien.py cp "$THIS_SCRIPT" alien://${MY_BINDIR}/${MY_JOBNAMEDATE}.sh@ALICE::CERN::EOS || true  # copy current job script to AliEn
#      [ ! "${CONTINUE_WORKDIR}" ] && alien.py cp "${MY_JOBSCRIPT}" alien://${MY_JOBWORKDIR}/alien_jobscript.sh@ALICE::CERN::EOS || true
    ) &> alienlog.txt

    pok "Submitting job \"${MY_JOBNAMEDATE}\" from $PWD"
    (
      echo "submit ${MY_JOBWORKDIR}/${MY_JOBNAMEDATE}.jdl" >> ${command_file}

      # finally we do a single call to alien:
      alien.py < ${command_file}
    ) &>> alienlog.txt

    MY_JOBID=$( (grep 'Your new job ID is' alienlog.txt | grep -oE '[0-9]+' || true) | sort -n | tail -n1)
    if [[ $MY_JOBID ]]; then
      pok "OK, display progress on https://alimonitor.cern.ch/agent/jobs/details.jsp?pid=$MY_JOBID"
    else
      per "Job submission failed: error log follows"
      cat alienlog.txt
    fi
  fi

  exit 0
fi

####################################################################################################
# The following part is executed on the worker node or locally
####################################################################################################
if [[ "${ONGRID}" == 0 ]]; then
  banner "Executing job in directory ${WORKDIR}"
  cd "${WORKDIR}" 2> /dev/null
fi

set -x


# ----------- START JOB PREAMBLE  ----------------------------- 
banner "Environment"
env

banner "OS detection"
lsb_release -a || true
cat /etc/os-release || true
cat /etc/redhat-release || true

if [ ! "$O2_ROOT" ]; then
  O2_PACKAGE_LATEST=`find /cvmfs/alice.cern.ch/el7-x86_64/Modules/modulefiles/O2 -type f -printf "%f\n" | tail -n1`
  banner "Loading O2 package $O2_PACKAGE_LATEST"
  [ "${O2TAG}" ] && O2_PACKAGE_LATEST=${O2TAG}
  eval "$(/cvmfs/alice.cern.ch/bin/alienv printenv O2::"$O2_PACKAGE_LATEST")"
fi
if [ ! "$XJALIEN_ROOT" ]; then
  XJALIEN_LATEST=`find /cvmfs/alice.cern.ch/el7-x86_64/Modules/modulefiles/xjalienfs -type f -printf "%f\n" | tail -n1`
  banner "Loading XJALIEN package $XJALIEN_LATEST"
  eval "$(/cvmfs/alice.cern.ch/bin/alienv printenv xjalienfs::"$XJALIEN_LATEST")"
fi
if [ ! "$O2DPG_ROOT" ]; then
  O2DPG_LATEST=`find /cvmfs/alice.cern.ch/el7-x86_64/Modules/modulefiles/O2DPG -type f -printf "%f\n" | tail -n1`
  banner "Loading O2DPG package $O2DPG_LATEST"
  eval "$(/cvmfs/alice.cern.ch/bin/alienv printenv O2DPG::"$O2DPG_LATEST")"
fi

banner "Running workflow"

# collect some common information
echo "CONT_WORKDIR ${CONTINUE_WORKDIR}"

cat /proc/cpuinfo > alien_cpuinfo.log 
cat /proc/meminfo > alien_meminfo.log

# ----------- PREPARE SOME ALIEN ENV -- useful for the job -----------

if [ "${ONGRID}" = "1" ]; then
  alien.py ps --jdl ${ALIEN_PROC_ID} > this_jdl.jdl
  ALIEN_JOB_OUTPUTDIR=$(grep "OutputDir" this_jdl.jdl | awk '//{print $3}' | sed 's/"//g' | sed 's/;//')
  export ALIEN_JOB_OUTPUTDIR

  # ----------- FETCH PREVIOUS CHECKPOINT IN CASE WE CONTINUE A JOB ----
  if [ "${CONTINUE_WORKDIR}" ]; then
    alien.py cp alien://${ALIEN_JOB_OUTPUTDIR}/checkpoint.tar.gz .
    tar -xzf checkpoint.tar.gz
    rm checkpoint.tar.gz
  fi
fi

# ----------- EXECUTE ACTUAL JOB  ------------------------------------ 
# source the actual job script from the work dir
chmod +x ./alien_jobscript.sh
./alien_jobscript.sh

# just to be sure that we get the logs
alien.py cp -f alien_log_${ALIEN_PROC_ID}.txt alien://${ALIEN_JOB_OUTPUTDIR}/lastlog.txt

# MOMENTARILU WE ZIP ALL LOG FILES
zip logs_PROCID${ALIEN_PROC_ID:-0}.zip *.log* alien_log*.txt

# We need to exit for the ALIEN JOB HANDLER!
exit 0
