#!/bin/bash

set -e

./scripts/clean-prefix.sh

{

./scripts/dist-update-hash.sh

./scripts/setenv.default.sh

./scripts/fftw.sh
./scripts/fftwf.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/qlat.sh

./scripts/c-lime.sh
./scripts/grid-sse4.sh
./scripts/gpt.sh

} |& tee $prefix/log.build.txt
