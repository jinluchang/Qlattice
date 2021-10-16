#!/bin/bash

set -e

./scripts/setenv.default.sh

./scripts/xz.sh
./scripts/tar.sh
./scripts/gmp.sh
./scripts/mpfr.sh
./scripts/mpc.sh
./scripts/gcc.sh
./scripts/perl.sh
./scripts/openssl.sh
./scripts/cmake.sh
./scripts/libffi.sh
./scripts/zlib.sh
./scripts/python.sh
./scripts/python-packages.sh
./scripts/re2c.sh
./scripts/ninja.sh
./scripts/llvm-project.sh

export CC=clang
export CXX=clang++

./scripts/openmpi.sh

export OMPI_CC=clang
export OMPI_CXX=clang++

./scripts/fftw.sh
./scripts/fftwf.sh
./scripts/cuba.sh
./scripts/eigen.sh
./scripts/qlat.sh

./scripts/autoconf.sh
./scripts/automake.sh
./scripts/c-lime.sh
./scripts/hdf5.sh
./scripts/grid.avx2.sh
./scripts/gpt.sh
./scripts/grid-tblum.avx2.sh
./scripts/hadrons-tblum.sh
