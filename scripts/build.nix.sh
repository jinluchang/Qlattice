#!/usr/bin/env bash

set -e

{
    ./scripts/setenv.default.sh
    ./scripts/qcore.sh

    ./scripts/qlat-utils.sh
    ./scripts/qlat.sh
    ./scripts/qlat-grid.sh
    ./scripts/qlat-cps.sh

    date
} 2>&1 | tee $prefix/log.build.txt
