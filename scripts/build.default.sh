#!/usr/bin/env bash

# `pipefail` so a failed step is not hidden by the `| tee` below (that masking
# is how a failed CPS build used to let this script exit 0 before qlat was
# installed).  CPS itself is optional for this configuration: it is not built
# here, and cps_python_package_optional() skips qlat-cps explicitly when the
# CPS library is absent (see scripts/cps-optional-functions.sh).
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

    ./scripts/qlat-utils.sh
    ./scripts/qlat.sh
    # The plain `default` flow does not build Grid (only `default-gpt` does), so
    # `grid-config` is legitimately absent here.  qlat-grid.sh fails loudly now
    # that it has `set -o pipefail`, so skip it explicitly instead of failing
    # the whole build (this is the same behaviour the old masked `| tee` gave).
    if which grid-config >/dev/null 2>&1 ; then
        ./scripts/qlat-grid.sh
    else
        echo "WARNING: Grid is not built in the 'default' configuration; skipping qlat-grid."
    fi
    cps_python_package_optional

    date
} 2>&1 | tee $prefix/log.build.txt
