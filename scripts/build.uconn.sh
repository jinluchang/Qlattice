#!/bin/bash

set -e

{

./scripts/setenv.uconn.sh
./scripts/qcore.sh

./scripts/gsl.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/perl.sh
./scripts/openssl.sh
./scripts/libffi.sh
./scripts/openblas.sh
./scripts/hdf5.sh
./scripts/python.sh
./scripts/python-pip.sh
./scripts/re2c.sh
./scripts/ninja.sh
./scripts/ninja-script.sh
./scripts/python-meson.sh
./scripts/python-meson-py.sh
./scripts/python-packages.sh

./scripts/fftw.sh

./scripts/qlat-utils.sh
./scripts/qlat.sh

./scripts/c-lime.sh
./scripts/gmp.sh
./scripts/mpfr.sh
./scripts/grid-clehner.avx2.sh
./scripts/qlat-grid.sh

./scripts/qmp.sh
./scripts/qio.sh
./scripts/cps.sh

./scripts/qlat-cps.sh

./scripts/gpt.sh

./scripts/qlat-examples-py.sh
./scripts/qlat-examples-cpp.sh
./scripts/qlat-examples-py-gpt.sh
./scripts/qlat-examples-py-cps.sh
./scripts/qlat-examples-cpp-grid.sh

} 2>&1 | tee $prefix/log.build.txt
