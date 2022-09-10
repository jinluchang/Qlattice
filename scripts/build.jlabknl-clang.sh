#!/bin/bash

set -e

{

./scripts/setenv.jlabknl.sh

export USE_COMPILER=gcc

./scripts/xz.sh
./scripts/tar.sh
./scripts/perl.sh
./scripts/openssl.sh
./scripts/libffi.sh
./scripts/zlib.sh
./scripts/openblas.sh
./scripts/python.sh
./scripts/python-packages.sh
./scripts/re2c.sh
./scripts/ninja.sh
./scripts/llvm-project.sh

export USE_COMPILER=clang

./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/autoconf.sh
./scripts/automake.sh
./scripts/c-lime.sh

./scripts/qlat.sh
./scripts/grid.knl.sh
./scripts/gpt.knl.sh

} |& tee $prefix/log.build.txt
