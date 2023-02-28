#!/bin/bash

export USE_COMPILER=clang

./scripts/qlat.clang.sh

name=qlat-grid

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    build="$prefix/build-clang"
    mkdir -p "$build"

    cd "$build"

    export CXX="$(grid-config --cxx)"
    # export CXX_LD="$(grid-config --cxxld)"

    touch "$wd"/qlat-grid/meson.build

    time-run meson setup "$wd/qlat-grid" \
        --prefix="$prefix"

    time-run meson compile -j$num_proc

    rm -rfv "$prefix"/lib
    time-run meson install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
