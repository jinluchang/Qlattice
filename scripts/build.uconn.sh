#!/bin/bash

# use command: fisbatch -n 4 -p generalsky

set -e

{

./scripts/setenv.uconn.sh

./scripts/gsl.sh
./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/perl.sh
./scripts/openssl.sh
./scripts/libffi.sh
./scripts/zlib.sh
./scripts/openblas.sh
./scripts/python.sh
./scripts/python-packages.sh
./scripts/qlat.sh

./scripts/c-lime.sh
./scripts/hdf5.sh
./scripts/grid.avx512.sh
./scripts/gpt.sh

} |& tee $prefix/log.build.txt
