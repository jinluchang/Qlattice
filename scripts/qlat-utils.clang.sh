#!/bin/bash

export USE_COMPILER=clang

. scripts/res/conf.sh

name=qlat-utils

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/build-qlat-utils-clang"
    mkdir -p "$build"

    cd "$build"

    rm -rfv "$prefix"/lib/python3*/*-packages/cqlat_utils.*
    rm -rfv "$prefix"/lib/python3*/*-packages/qlat_utils

    if [ -n "$QLAT_CXX" ] ; then
        export CXX="$QLAT_CXX"
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

    touch "$wd"/qlat-utils/meson.build

    prefix_python="$prefix/lib/python3/qlat-packages"

    meson setup "$wd/qlat-utils" \
        -Dpython.platlibdir="$prefix_python" -Dpython.purelibdir="$prefix_python" \
        --prefix="$prefix"
    time meson compile -j$num_proc
    meson install

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
