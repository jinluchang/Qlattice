#!/bin/bash

set -e

{

./scripts/setenv.default.sh
./scripts/qcore.sh

./scripts/ninja-script.sh

./scripts/qlat-utils.sh
./scripts/qlat.sh
./scripts/qlat-grid.sh

} 2>&1 | tee $prefix/log.build.txt
