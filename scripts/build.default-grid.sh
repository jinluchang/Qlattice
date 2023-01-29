#!/bin/bash

set -e

{

./scripts/setenv.default.sh

./scripts/c-lime.sh
./scripts/grid.avx2.sh

} 2>&1 | tee $prefix/log.build.txt
