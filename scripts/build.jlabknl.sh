#!/bin/bash

set -e

{

./scripts/setenv.jlabknl.sh

export CC=gcc
export CXX=g++
export CFLAGS=" "
export CXXFLAGS=" "

./scripts/xz.sh
./scripts/tar.sh
./scripts/perl.sh
./scripts/openssl.sh
./scripts/libffi.sh
./scripts/zlib.sh
./scripts/openblas.sh
./scripts/python.sh
./scripts/python-packages.sh

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
export CFLAGS=
export CXXFLAGS=

./scripts/qlat.sh
./scripts/grid.knl.sh
./scripts/gpt.knl.sh

} |& tee $prefix/log.build.txt
