#!/bin/bash

./scripts/qlat.sh

. scripts/conf.sh

name=qlat-grid

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/build-qlat-grid"
    mkdir -p "$build"

    cd "$build"

    rm -rfv "$prefix"/lib/python3*/*-packages/cqlat_grid.*
    rm -rfv "$prefix"/lib/python3*/*-packages/qlat_grid

    export CXX="$(grid-config --cxx)"
    # export CXX_LD="$(grid-config --cxxld)"

    ( cd "$wd/qlat-grid/pylib/cqlat_grid" ; bash update.sh )

    touch "$wd"/qlat-grid/meson.build

    prefix_python="$prefix/lib/python3/qlat-packages"

    meson setup "$wd/qlat-grid" \
        -Dpython.platlibdir="$prefix_python" -Dpython.purelibdir="$prefix_python" \
        --prefix="$prefix"
    ninja -j$num_proc
    ninja install

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
