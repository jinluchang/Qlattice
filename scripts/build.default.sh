#!/bin/bash

set -e

{

./scripts/setenv.default.sh
./scripts/qcore.sh

./scripts/cuba.sh
./scripts/ninja-script.sh

./scripts/qlat-grid.sh

} 2>&1 | tee $prefix/log.build.txt
