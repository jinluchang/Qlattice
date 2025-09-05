#!/usr/bin/env bash

name=setenv-perlmutter-cpu

source qcore/set-prefix.sh

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh

#
cat >"$prefix/setenv.sh" <<EOF
if [ -z "\$num_proc" ] ; then
    export num_proc=16
fi

module load PrgEnv-gnu


#module purge

module load craype-x86-milan

# System-provided math & I/O libraries
module load cray-libsci           # BLAS/LAPACK/ScaLAPACK
module load cray-fftw             # FFTW (parallel)
module load cray-hdf5-parallel    # MPI-enabled HDF5

# Python toolchain
module load python                # Cray-packaged Python + mpi4py

export PYTHONPATH="\$prefix/python-packages/lib/python3.11/site-packages:\$PYTHONPATH"

# Support modules
module load xpmem                 # XPMEM for shared-memory transport
module load perftools-base        # Cray performance tools


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
export QLAT_CXX=\${GCC_PATH}/g++
export QLAT_MPICXX=mpic++
EOF

    #
    "$wd"/qcore/bin/mk-setenv-dir.sh --keep
    echo "!!!! $name build !!!!"
} } 2>&1 | tee $prefix/log.$name.txt
