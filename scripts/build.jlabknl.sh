#!/bin/bash

set -e

{

./scripts/setenv.jlabknl.sh

# export USE_COMPILER=gcc

./scripts/xz.sh
./scripts/tar.sh
./scripts/perl.sh
./scripts/openssl.sh
./scripts/libffi.sh
./scripts/zlib.sh
./scripts/openblas.sh
./scripts/python.sh
./scripts/python-packages.sh

./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/autoconf.sh
./scripts/automake.sh
./scripts/c-lime.sh
./scripts/hdf5.sh

# ./scripts/openmpi.sh

./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh

# export USE_COMPILER=intel

./scripts/qlat.sh
./scripts/grid.knl.sh
./scripts/gpt.sh

} |& tee $prefix/log.build.txt
