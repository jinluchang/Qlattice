#!/bin/bash

name=qlat-docs

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p "$build_dir"
    rsync -a --delete "$wd"/docs/ "$build_dir"/

    cd "$build_dir"

    rm -rfv generated _build

    make html latexpdf

    mkdir -p "$prefix"/share/doc/qlat

    rsync -a --delete "$build_dir"/_build/ "$prefix"/share/doc/qlat

    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
