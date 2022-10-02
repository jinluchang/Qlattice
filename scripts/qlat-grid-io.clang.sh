#!/bin/bash

./scripts/qlat.clang.sh

export USE_COMPILER=clang

. scripts/conf.sh

name=qlat-grid-io

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/build-qlat-grid-io-clang"
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

} |& tee $prefix/log.$name.txt
