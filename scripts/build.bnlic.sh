#!/bin/bash

set -e

{

./scripts/setenv.bnlic.sh

./scripts/tar.sh

./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/ninja.sh
./scripts/python-meson.sh

./scripts/qlat.sh

./scripts/c-lime.sh
./scripts/hdf5.sh
./scripts/grid.avx2.sh
./scripts/gpt.sh
./scripts/qlat-grid-io.sh

} |& tee $prefix/log.build.txt
