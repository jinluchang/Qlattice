#!/bin/bash

set -e

{

./scripts/setenv.bnlic.sh
./scripts/qcore.sh

./scripts/tar.sh
./scripts/zlib.sh
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

./scripts/fftw_mpi.sh
./scripts/cuba.sh
./scripts/eigen.sh

./scripts/qlat.sh

./scripts/qlat-examples-py.sh
./scripts/qlat-examples-cpp.sh

./scripts/c-lime.sh
./scripts/hdf5.sh
./scripts/grid-clehner.avx2.sh
./scripts/gpt.sh
./scripts/qlat-grid.sh

./scripts/qlat-examples-py-gpt.sh

} 2>&1 | tee $prefix/log.build.txt
