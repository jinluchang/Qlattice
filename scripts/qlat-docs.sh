#!/bin/bash

name=qlat-docs

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p "$build_dir"
    time-run rsync -a --delete "$wd"/docs/ "$build_dir"/

    cd "$build_dir"

    rm -rfv generated _build

    time-run make html || true
    time-run make latexpdf || true

    mkdir -p "$prefix"/share/doc/qlat

    time-run rsync -a --delete "$build_dir"/_build/ "$prefix"/share/doc/qlat

    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
