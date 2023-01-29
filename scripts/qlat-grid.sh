#!/bin/bash

./scripts/qlat.sh

name=qlat-grid

source qcore/set-prefix.sh $name

{ time {

    echo "!!!! build $name !!!!"

    build="$prefix/build"
    mkdir -p "$build"

    cd "$build"

    rm -rfv "$prefix"/lib/python3

    export CXX="$(grid-config --cxx)"
    # export CXX_LD="$(grid-config --cxxld)"

    touch "$wd"/qlat-grid/meson.build

    prefix_python="$prefix/lib/python3"

    meson setup "$wd/qlat-grid" \
        -Dpython.platlibdir="$prefix_python" -Dpython.purelibdir="$prefix_python" \
        --prefix="$prefix"

    time meson compile -j$num_proc
    meson install

    cd "$wd"

    mk-setenv.sh

    echo "!!!! $name build !!!!"

    rm -rf "$temp_dir" || true

} } 2>&1 | tee $prefix/log.$name.txt
