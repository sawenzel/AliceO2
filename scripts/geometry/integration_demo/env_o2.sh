# Environment for everything that links O2 (o2-sim, root, the harness).
# The eval and every command using it must run in the SAME shell invocation.
export ALIBUILD_WORK_DIR=$HOME/alisw/sw
export B=$HOME/alisw/sw/BUILD/O2-latest-swenzel-bvhsurfacesolid/O2
export SW=$HOME/alisw/sw/ubuntu2404_aarch64
cd $B && eval "$($HOME/alisw/alibuild/alienv printenv O2/latest-swenzel-bvhsurfacesolid-o2,ninja/latest,CMake/latest)"
export LD_LIBRARY_PATH=$B/stage/lib:$B/stage/lib64:$LD_LIBRARY_PATH
export VMCWORKDIR=$HOME/alisw/O2
export ROOT_INCLUDE_PATH=$HOME/alisw/O2/Detectors/Base/include:$ROOT_INCLUDE_PATH
