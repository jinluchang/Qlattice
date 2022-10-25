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

    if [ -n "$QLAT_MPICXX" ] ; then
        export CXX="$QLAT_MPICXX"
        export MPICXX=false
    fi
    if [ -n "$QLAT_CXXFLAGS" ] ; then
        export CXXFLAGS="$QLAT_CXXFLAGS"
    fi
    if [ -n "$QLAT_LDFLAGS" ] ; then
        export LDFLAGS="$QLAT_LDFLAGS"
    fi
    if [ -n "$QLAT_LIBS" ] ; then
        export LIBS="$QLAT_LIBS"
    fi

    touch "$wd"/qlat-grid-io/meson.build

    meson "$wd/qlat-grid-io" --prefix="$prefix" -Dgrid_prefix="$prefix"
    ninja -j$num_proc
    ninja install

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
