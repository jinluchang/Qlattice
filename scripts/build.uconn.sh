#!/bin/bash

# use command: fisbatch -n 4 -p generalsky

set -e

{

./scripts/setenv.uconn.sh

./scripts/gsl.sh
./scripts/fftw.sh
./scripts/fftwf.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/qlat.sh

./scripts/c-lime.sh
./scripts/grid.avx512.sh
./scripts/gpt.sh

} |& tee $prefix/log.build.txt
