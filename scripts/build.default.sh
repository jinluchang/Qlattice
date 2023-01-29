#!/bin/bash

set -e

{

./scripts/setenv.default.sh

./scripts/qcore.sh

./scripts/ninja-script.sh

./scripts/qlat.sh

} 2>&1 | tee $prefix/log.build.txt
