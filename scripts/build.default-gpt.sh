#!/bin/bash

set -e

{

./scripts/setenv.default.sh
./scripts/qcore.sh

./scripts/cuba.sh
./scripts/eigen.sh
./scripts/ninja-script.sh

./scripts/qlat-utils.sh
./scripts/qlat.sh

./scripts/c-lime.sh
./scripts/grid-clehner.avx2.sh
./scripts/gpt.sh
./scripts/qlat-grid.sh

} 2>&1 | tee $prefix/log.build.txt
