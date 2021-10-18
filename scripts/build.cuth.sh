#!/bin/bash

set -e

./scripts/clean-prefix.sh

{

./scripts/dist-update-hash.sh

./scripts/setenv.default.sh

./scripts/xz.sh
./scripts/tar.sh
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
./scripts/python.sh
./scripts/python-packages.sh
./scripts/re2c.sh
./scripts/ninja.sh
./scripts/llvm-project.sh

./scripts/openmpi.sh

./scripts/fftw.sh
./scripts/fftwf.sh
./scripts/cuba.sh
./scripts/eigen.sh
./scripts/qlat.sh

./scripts/autoconf.sh
./scripts/automake.sh
./scripts/c-lime.sh
./scripts/hdf5.sh
./scripts/grid.gen16.sh
./scripts/gpt.sh

} |& tee $prefix/log.build.txt
