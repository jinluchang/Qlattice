#!/bin/bash

# Usage:
# find-library.sh lib-name-1 lib-name-2 ...
# return first prefix of the path in LIBRARY_PATH (endswith /lib or /lib64) that contains one of the lib-names

if python-check-version.py >/dev/null 2>&1 && find-library.py "$@" ; then
    # successfully finished.
    exit 0
fi

# TODO

echo ""
exit 1
