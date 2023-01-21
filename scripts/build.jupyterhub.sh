#!/bin/bash

# Ubuntu packages: libgsl-dev zlib1g-dev libssl-dev libeigen3-dev libopenmpi-dev libfftw3-dev libfftw3-mpi-dev libhdf5-dev libhdf5-mpi-dev  libmpfr-dev ninja-build libsqlite3-dev libgdbm-dev liblzma-dev libbz2-dev uuid-dev tk-dev libgeos-dev texlive-metapost

set -e

{

./scripts/setenv.default.sh

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

./scripts/qlat-examples-py.sh
./scripts/qlat-examples-cpp.sh

./scripts/c-lime.sh
./scripts/grid.avx2.sh
./scripts/gpt.sh

./scripts/qlat-packages.sh

./scripts/qlat-examples-py-gpt.sh

} 2>&1 | tee $prefix/log.build.txt
