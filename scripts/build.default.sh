#!/bin/bash

set -e

{

./scripts/setenv.default.sh
./scripts/qcore.sh

./scripts/ninja-script.sh
./scripts/eigen.sh

./scripts/qlat-utils.sh
./scripts/qlat.sh
./scripts/qlat-grid.sh
./scripts/qlat-cps.sh

} 2>&1 | tee $prefix/log.build.txt
