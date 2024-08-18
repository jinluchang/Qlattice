#!/usr/bin/env bash

./scripts/update-sources.sh

name=dist-update-hash

source qcore/set-prefix.sh $name

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh

    distfiles="$wd/distfiles"

    cd "$distfiles"

    cat sha256sums.txt

} } 2>&1 | tee $prefix/log.$name.txt

