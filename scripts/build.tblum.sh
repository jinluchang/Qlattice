#!/bin/bash

set -e

./scripts/setenv.sh

./scripts/fftw.sh
./scripts/fftwf.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/qlat.sh

./scripts/c-lime.sh
./scripts/hdf5.sh
./scripts/grid-tblum-avx2.sh
./scripts/hadrons-tblum.sh
./scripts/grid-avx2.sh
./scripts/gpt.sh
