# Environment for O2_CADtoTGeo.py and any other pythonOCC script.
# NEVER source this into a shell that runs o2-sim: the pythonOCC PYTHONPATH
# prepends make o2-sim segfault at startup.
SW=$HOME/alisw/sw/ubuntu2404_aarch64
export LD_LIBRARY_PATH=$SW/OCCT/latest/lib:$SW/Python/latest/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$SW/pythonOCC/latest/lib/python3.10/site-packages:$SW/Python-modules/latest/lib/python3.10/site-packages:$PYTHONPATH
PYOCC=$SW/Python/latest/bin/python3.10
