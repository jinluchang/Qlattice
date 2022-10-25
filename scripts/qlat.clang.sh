#!/bin/bash

./scripts/qlat-utils.clang.sh

export USE_COMPILER=clang

. scripts/conf.sh

name=qlat

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/build-qlat-clang"
    mkdir -p "$build"

    cd "$build"

    rm -rfv "$prefix"/include/qlat
    rm -rfv "$prefix"/include/qlat-setup.h
    rm -rfv "$prefix"/lib/python3*/*-packages/cqlat.*
    rm -rfv "$prefix"/lib/python3*/*-packages/qlat
    rm -rfv "$prefix"/lib/python3*/*-packages/qlat_gpt.py
    rm -rfv "$prefix"/lib/python3*/*-packages/rbc_ukqcd*
    rm -rfv "$prefix"/lib/python3*/*-packages/auto_contractor*

    if [ -n "$QLAT_MPICXX" ] ; then
        export CXX="$QLAT_MPICXX"
        export MPICXX=false
        option="-Duse_cxx=true"
    else
        option=
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

    touch "$wd"/qlat/meson.build

    prefix_python="$prefix/lib/python3/qlat-packages"

    meson "$wd/qlat" \
        -Dpython.platlibdir="$prefix_python" -Dpython.purelibdir="$prefix_python" \
        --prefix="$prefix" $option
    ninja -j$num_proc
    ninja install

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
