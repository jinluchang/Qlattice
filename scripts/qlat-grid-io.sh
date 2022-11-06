#!/bin/bash

./scripts/qlat.sh

. scripts/conf.sh

name=qlat-grid-io

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/build-qlat-grid-io"
    mkdir -p "$build"

    cd "$build"

    rm -rfv "$prefix"/include/qlat-grid-io
    rm -rfv "$prefix"/lib/python3*/*-packages/cqlat-grid-io.*
    rm -rfv "$prefix"/lib/python3*/*-packages/qlat-grid-io

    export CXX="$(grid-config --cxx)"
    # export CXX_LD="$(grid-config --cxxld)"

    ( cd "$wd/qlat-grid-io/pylib/cqlat_grid_io" ; bash update.sh )

    touch "$wd"/qlat-grid-io/meson.build

    prefix_python="$prefix/lib/python3/qlat-packages"

    meson "$wd/qlat-grid-io" \
        -Dpython.platlibdir="$prefix_python" -Dpython.purelibdir="$prefix_python" \
        --prefix="$prefix"
    ninja -j$num_proc
    ninja install

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
