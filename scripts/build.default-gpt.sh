#!/bin/bash

# Need to have gcc, mpi, python3 (numpy, simpy, psutil, meson), openssl, ninja, zlib, gsl installed beforehand.

# Ubuntu packages: libgsl-dev zlib1g-dev libssl-dev libopenmpi-dev python3-sympy python3-numpy python3-scipy python3-psutil libmpfr-dev meson ninja-build

set -e

{

./scripts/setenv.default.sh

./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh
./scripts/cuba.sh
./scripts/eigen.sh

./scripts/qlat.sh

./scripts/c-lime.sh
./scripts/grid.avx2.sh
./scripts/gpt.sh
./scripts/qlat-grid-io.sh

} 2>&1 | tee $prefix/log.build.txt
