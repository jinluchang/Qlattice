#!/bin/bash

set -e

# need to have mpi, python3 (numpy, simpy), openssl installed beforehand.

./scripts/setenv.default.sh

./scripts/fftw.sh
./scripts/fftwf.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/qlat.sh

./scripts/c-lime.sh
./scripts/grid.avx2.sh
./scripts/gpt.sh
