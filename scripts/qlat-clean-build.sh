#!/usr/bin/env bash

name=qlat-clean-build

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh

    rm -rfv "$prefix/../qlat-utils/build"
    rm -rfv "$prefix/../qlat/build"
    rm -rfv "$prefix/../qlat-grid/build"
    rm -rfv "$prefix/../qlat-cps/build"

    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
