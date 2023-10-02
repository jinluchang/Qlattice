#!/bin/bash

set -e

{
    ./scripts/setenv.frontier.sh

    ./scripts/grid-clehner.frontier.sh
    ./scripts/gpt.sh

    ./scripts/qlat-utils.sh
    ./scripts/qlat.sh
    ./scripts/qlat-grid.sh
    ./scripts/qlat-cps.sh

    date
} 2>&1 | tee $prefix/log.build.txt
