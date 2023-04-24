#!/bin/bash

# Ubuntu packages: libgsl-dev zlib1g-dev libssl-dev libeigen3-dev libopenmpi-dev libfftw3-dev libfftw3-mpi-dev libmpfr-dev ninja-build libsqlite3-dev libgdbm-dev liblzma-dev libbz2-dev uuid-dev tk-dev libgeos-dev texlive-metapost libopenblas-dev

set -e

{

rm -rfv ~/.cache/pip

./scripts/setenv.default.sh
./scripts/qcore.sh

./scripts/hdf5.sh

./scripts/python.sh
./scripts/python-pip.sh
./scripts/ninja-script.sh
./scripts/python-meson.sh
./scripts/python-meson-py.sh
./scripts/python-packages.sh
./scripts/python-jupyter.sh

./scripts/cuba.sh
./scripts/eigen.sh

./scripts/qlat-packages.sh

./scripts/c-lime.sh
./scripts/qmp.sh
./scripts/qio.sh
./scripts/cps.sh
./scripts/grid-clehner.avx2.sh
./scripts/gpt.sh

./scripts/qlat-packages.sh

./scripts/qlat-examples-py.sh
./scripts/qlat-examples-cpp.sh
./scripts/qlat-examples-py-gpt.sh
./scripts/qlat-examples-py-cps.sh
./scripts/qlat-examples-cpp-grid.sh

} 2>&1 | tee $prefix/log.build.txt
