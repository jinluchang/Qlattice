#!/bin/bash

# Need to have gcc, eigen3, mpi, fftw3-mpi, python3 (numpy, simpy, psutil, meson), openssl, ninja, zlib, gsl installed beforehand.

# Ubuntu packages: libgsl-dev zlib1g-dev libssl-dev libeigen3-dev libopenmpi-dev libfftw3-mpi-dev python3-sympy python3-numpy python3-scipy python3-psutil libmpfr-dev meson ninja-build

set -e

{

./scripts/setenv.default.sh

./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh
./scripts/cuba.sh
./scripts/eigen.sh

./scripts/qlat.sh

} 2>&1 | tee $prefix/log.build.txt
