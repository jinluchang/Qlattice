#!/bin/bash

set -e

./scripts/clean-prefix.sh

{

./scripts/dist-update-hash.sh

./scripts/setenv.jlabknl.sh

export CC=gcc
export CXX=g++

./scripts/xz.sh
./scripts/tar.sh
./scripts/perl.sh
./scripts/openssl.sh
./scripts/libffi.sh
./scripts/zlib.sh
./scripts/python.sh
./scripts/python-packages.sh

./scripts/fftw.sh
./scripts/fftwf.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/autoconf.sh
./scripts/automake.sh
./scripts/c-lime.sh
./scripts/hdf5.sh

export CC=
export CXX=

./scripts/qlat.sh
./scripts/grid-tblum.knl.sh
./scripts/hadrons-tblum.sh
./scripts/grid.knl.sh
./scripts/gpt.sh

} |& tee $prefix/log.build.txt
