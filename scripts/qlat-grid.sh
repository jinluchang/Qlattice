#!/bin/bash

./scripts/qlat.sh

name=qlat-grid

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    build="$prefix/build"
    mkdir -p "$build"

    cd "$build"

    rm -rfv "$prefix"/lib/python3

    export CXX="$(grid-config --cxx)"
    # export CXX_LD="$(grid-config --cxxld)"

    touch "$wd"/qlat-grid/meson.build

    meson setup "$wd/qlat-grid" \
        --prefix="$prefix"

    time meson compile -j$num_proc
    meson install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
