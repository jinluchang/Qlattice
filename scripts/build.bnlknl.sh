#!/bin/bash

set -e

{

./scripts/setenv.bnlknl.sh

export CC=gcc
export CXX=g++

./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/autoconf.sh
./scripts/automake.sh
./scripts/c-lime.sh

export CC=
export CXX=

./scripts/qlat.sh
./scripts/grid.knl.sh
./scripts/gpt.knl.sh

} |& tee $prefix/log.build.txt
