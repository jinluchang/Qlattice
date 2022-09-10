#!/bin/bash

set -e

{

./scripts/setenv.bnlknl.sh

export USE_COMPILER=gcc

./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/autoconf.sh
./scripts/automake.sh
./scripts/c-lime.sh

export USE_COMPILER=intel

./scripts/qlat.sh
./scripts/grid.knl.sh
./scripts/gpt.knl.sh

} |& tee $prefix/log.build.txt
