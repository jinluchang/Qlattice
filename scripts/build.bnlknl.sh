#!/bin/bash

set -e

{

./scripts/setenv.bnlknl.sh

export CC=gcc
export CXX=g++

./scripts/fftw.sh
./scripts/fftwf.sh
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
./scripts/gpt.sh

} |& tee $prefix/log.build.txt
