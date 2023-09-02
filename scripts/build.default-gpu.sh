#!/bin/bash

set -e

{
    ./scripts/setenv.gpu.sh
    ./scripts/qcore.sh

    ./scripts/python-venv.sh
    ./scripts/python-pip-install.sh

    ./scripts/cuba.sh
    ./scripts/eigen.sh
    ./scripts/ninja-script.sh

    ./scripts/c-lime.sh
    ./scripts/qmp.sh
    ./scripts/qio.sh
    ./scripts/cps.sh
    ./scripts/grid-clehner.gpu.sh
    ./scripts/gpt.sh

    ./scripts/qlat-utils.sh
    ./scripts/qlat.sh
    ./scripts/qlat-grid.sh
    ./scripts/qlat-cps.sh

    date
} 2>&1 | tee $prefix/log.build.txt
