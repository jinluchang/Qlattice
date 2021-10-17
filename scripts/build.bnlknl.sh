#!/bin/bash

set -e

./scripts/setenv.bnlknl.sh

./scripts/fftw.sh
./scripts/fftwf.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/autoconf.sh
./scripts/automake.sh
./scripts/c-lime.sh

./scripts/qlat.sh
./scripts/grid.knl.sh
./scripts/gpt.sh
