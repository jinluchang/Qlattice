#!/bin/bash

# need to have gcc, mpi, python3 (numpy, simpy, psutil, meson), openssl, ninja installed beforehand.

set -e

{

./scripts/setenv.default.sh

./scripts/gsl.sh
./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/qlat-utils.sh
./scripts/qlat.sh

} |& tee $prefix/log.build.txt
