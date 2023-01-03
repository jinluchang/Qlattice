#!/bin/bash

set -e

{

./scripts/setenv.bnlknl.sh

./scripts/tar.sh

./scripts/zlib.sh
./scripts/xz.sh
./scripts/gsl.sh
./scripts/gmp.sh
./scripts/mpfr.sh
./scripts/mpc.sh
./scripts/isl.sh
./scripts/gcc.sh
./scripts/binutils.sh
./scripts/perl.sh
./scripts/openssl.sh
./scripts/libffi.sh
./scripts/openblas.sh
./scripts/python.sh
./scripts/python-pip.sh
./scripts/re2c.sh
./scripts/ninja.sh
./scripts/ninja-script.sh
./scripts/python-meson.sh
./scripts/python-meson-py.sh
./scripts/python-packages.sh

./scripts/cuba.sh
./scripts/eigen.sh

./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh

./scripts/qlat.sh

./scripts/qlat-examples-py.sh
./scripts/qlat-examples-cpp.sh

./scripts/autoconf.sh
./scripts/automake.sh
./scripts/c-lime.sh
./scripts/hdf5.sh
./scripts/grid.knl.sh
./scripts/gpt.sh
./scripts/qlat-grid.sh

./scripts/qlat-examples-py-gpt.sh

} 2>&1 | tee $prefix/log.build.txt
