#!/usr/bin/env bash

name=setenv-frontier

# https://github.com/paboyle/Grid/blob/develop/systems/Frontier/benchmarks/bench2.slurm

source qcore/set-prefix.sh

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh

#
cat >"$prefix/setenv.sh" <<EOF
if [ -z "\$num_proc" ] ; then
    export num_proc=8
fi
# module purge
#
if [ "\$PRG_ENV" = GNU ] ; then
    module load PrgEnv-gnu
else
    module load PrgEnv-amd
fi
#
module load craype-accel-amd-gfx90a
module load rocm
#
export MPICH_GPU_SUPPORT_ENABLED=1
export CFLAGS="-I\${ROCM_PATH}/include"
export CXXFLAGS="-I\${ROCM_PATH}/include"
export LDFLAGS="-L\${ROCM_PATH}/lib -lamdhip64"
#
module load craype-x86-trento
module load xpmem
module load perftools-base
module load openblas
#
module list
#
if [ -z "\$USE_COMPILER" ] ; then
    export USE_COMPILER=gcc
fi
#
export CC=cc
export CXX=CC
export MPICC=cc
export MPICXX=CC
export QLAT_CXX=amdclang++
export QLAT_MPICXX=mpic++
EOF

    #
    "$wd"/qcore/bin/mk-setenv-dir.sh --keep
    echo "!!!! $name build !!!!"
} } 2>&1 | tee $prefix/log.$name.txt
