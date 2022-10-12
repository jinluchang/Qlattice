#!/bin/bash

set -e

{

./scripts/setenv.jlabknl.sh

./scripts/xz.sh
./scripts/tar.sh
./scripts/perl.sh
./scripts/openssl.sh
./scripts/zlib.sh
./scripts/fftw_mpi.sh
./scripts/fftwf_mpi.sh

./scripts/autoconf.sh
./scripts/automake.sh
./scripts/c-lime.sh
./scripts/hdf5.sh

export USE_COMPILER="intel"

./scripts/grid-tblum.knl.sh
./scripts/hadrons-tblum.sh

export USE_COMPILER="gcc"

./scripts/cuba.sh
./scripts/eigen.sh

./scripts/libffi.sh
./scripts/openblas.sh
./scripts/python.sh
./scripts/python-pip.sh
./scripts/ninja.sh
./scripts/python-meson.sh
./scripts/python-packages.sh

export USE_COMPILER="intel"

./scripts/qlat.sh
./scripts/grid.knl.sh
./scripts/gpt.sh
./scripts/qlat-grid-io.sh

} |& tee $prefix/log.build.txt
