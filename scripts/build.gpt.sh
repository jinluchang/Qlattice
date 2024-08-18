#!/usr/bin/env bash

set -e

{
    ./scripts/setenv.default.sh
    ./scripts/qcore.sh

    ./scripts/cuba.sh
    ./scripts/eigen.sh
    ./scripts/ninja-script.sh

    ./scripts/c-lime.sh
    ./scripts/qmp.sh
    ./scripts/qio.sh
    ./scripts/cps.sh
    if cat /proc/cpuinfo | grep avx2 ; then
        ./scripts/grid-clehner.avx2.sh
    else
        ./scripts/grid-clehner.gen16.sh
    fi
    ./scripts/gpt.sh

    date
} 2>&1 | tee $prefix/log.build.txt
