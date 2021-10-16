#!/bin/bash

set -e

./scripts/setenv.jlabknl.sh

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

export CXX=icpc
export MPICXX=mpiicpc

./scripts/qlat.sh
./scripts/grid-tblum-knl.sh
./scripts/hadrons-tblum.sh
./scripts/grid-knl.sh
./scripts/gpt.sh
