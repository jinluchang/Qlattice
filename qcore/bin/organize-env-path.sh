#!/bin/bash

# Usage:
# if which organize-env-path.sh >/dev/null 2>&1 ; then
#     source <(organize-env-path.sh)
# fi

if ! python-check-version.py >/dev/null 2>&1 ; then
    echo "echo 'python3 does not work for organize-env-path.py. Don't do anything.'"
    exit 1
fi

organize-env-path.py "$@"
