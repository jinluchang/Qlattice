#!/bin/bash

set -e

{

./scripts/setenv.ofpknl.sh

# export CC=gcc
# export CXX=g++

./scripts/zlib.sh
# ./scripts/perl.sh
# ./scripts/openssl.sh
./scripts/libffi.sh
# ./scripts/openblas.sh
./scripts/python.sh
./scripts/python-packages.sh

./scripts/fftw.sh
./scripts/fftwf.sh
./scripts/cuba.sh
./scripts/eigen.sh
./scripts/c-lime.sh

# export CC=
# export CXX=

./scripts/qlat.sh
./scripts/grid.knl.sh
./scripts/gpt.sh

} |& tee $prefix/log.build.txt
