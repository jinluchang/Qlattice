#!/usr/bin/env bash

set -e

{

    ./scripts/setenv.default.sh

} 2>&1 | tee $prefix/log.build.txt
