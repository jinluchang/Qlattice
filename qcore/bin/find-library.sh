#!/usr/bin/env bash

# Usage:
# find-library.sh lib-name-1 lib-name-2 ...
# return first prefix of the path in LIBRARY_PATH (endswith /lib or /lib64) that contains one of the lib-names

if python-check-version.py >/dev/null 2>&1 && find-library.py "$@" ; then
    # successfully finished.
    exit 0
fi

# Fall back to the bash version (only deal with 1 argument case)

libname="$1"

value="$LIBRARY_PATH"
if [ -n "$value" ] ; then
    IFS=':' read -a vs <<< "$value"
    for v in "${vs[@]}" ; do
        if [ -f "$v"/"$libname"* ] ; then
            echo "$(dirname "$v")"
            exit 0
        fi
    done
fi
