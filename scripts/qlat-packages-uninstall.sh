#!/usr/bin/env bash

name=qlat-packages-uninstall

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    time-run pip3 uninstall -y qlat-utils qlat qlat-grid qlat-cps || true

    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
