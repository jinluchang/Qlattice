#!/usr/bin/env bash

# `pipefail` so a failed step is not hidden by the `| tee` below (that masking
# is how a failed CPS build used to let this script exit 0 before qlat was
# installed).  CPS is legitimately unavailable on macOS (missing GSL headers),
# so the CPS build is wrapped as optional there; on Linux a failed CPS build is
# still fatal (see scripts/cps-optional-functions.sh).
set -e -o pipefail

if ! [ -f ./scripts/cps-optional-functions.sh ] ; then
    echo "ERROR: ./scripts/cps-optional-functions.sh not found; run this script from the repository root." >&2
    exit 1
fi
source ./scripts/cps-optional-functions.sh

{
    ./scripts/setenv.default.sh
    ./scripts/qcore.sh

    ./scripts/python-venv.sh
    ./scripts/python-pip-install.sh

    ./scripts/cuba.sh
    ./scripts/eigen.sh
    ./scripts/ninja-script.sh

    ./scripts/c-lime.sh
    ./scripts/qmp.sh
    ./scripts/qio.sh
    cps_build_optional cps ./scripts/cps.sh
    if cat /proc/cpuinfo | grep avx2 ; then
        ./scripts/grid-clehner.avx2.sh
    else
        ./scripts/grid-clehner.gen16.sh
    fi
    ./scripts/gpt.sh

    ./scripts/qlat-utils.sh
    ./scripts/qlat.sh
    ./scripts/qlat-grid.sh
    cps_python_package_optional

    date
} 2>&1 | tee $prefix/log.build.txt
