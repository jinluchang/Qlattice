#!/bin/bash

set -e

{

./scripts/setenv.default.sh

./scripts/xz.sh
./scripts/tar.sh
./scripts/gsl.sh
./scripts/gmp.sh
./scripts/mpfr.sh
./scripts/mpc.sh
./scripts/gcc.sh
./scripts/binutils.sh
./scripts/perl.sh
./scripts/openssl.sh
./scripts/cmake.sh
./scripts/libffi.sh
./scripts/zlib.sh
./scripts/openblas.sh
./scripts/python.sh
./scripts/python-packages.sh
./scripts/re2c.sh
./scripts/ninja.sh
./scripts/llvm-project.sh

./scripts/openmpi.sh

./scripts/fftw.sh
./scripts/fftwf.sh
./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh
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

} |& tee $prefix/log.build.txt
