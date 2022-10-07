#!/bin/bash

# Need to run on computing node for UCONN HPC
# You can use within interactive session obtained by the following command:
# fisbatch -n 2 -p generalsky

set -e

{

./scripts/setenv.uconn.sh

./scripts/gsl.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/perl.sh
./scripts/openssl.sh
./scripts/libffi.sh
./scripts/openblas.sh
./scripts/python.sh
./scripts/python-pip.sh
./scripts/ninja.sh
./scripts/python-meson.sh
./scripts/python-packages.sh

./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh

./scripts/qlat.sh

./scripts/c-lime.sh
./scripts/hdf5.sh
./scripts/grid.avx512.sh
./scripts/gpt.sh
./scripts/qlat-grid-io.sh

} |& tee $prefix/log.build.txt
