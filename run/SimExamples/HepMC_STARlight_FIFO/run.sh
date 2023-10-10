#!/usr/bin/env bash
#
# This is a simple simulation example showing the following things
#

# a) how to start STARlight event generator in a clean environment to produce HepMC files
# b) how to run simulation reading from the HepMC file

set -x

# PART a)
#env -i HOME="$HOME" USER="$USER" PATH="/bin:/usr/bin:/usr/local/bin" \
#    ALIBUILD_WORK_DIR="$ALIBUILD_WORK_DIR" 

# PART a)
# ./run-starlight.sh

# PART b)
NEV=10
o2-sim -j 4 -n ${NEV} -g hepmc -m PIPE ITS -o sim \
       --configKeyValues "FileOrCmd.cmd=${PWD}/starlight.sh;Diamond.position[2]=0.1;Diamond.width[2]=0.05"
