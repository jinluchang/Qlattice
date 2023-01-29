#!/bin/bash

./scripts/qlat-utils.sh

name=qlat

source qcore/set-prefix.sh $name

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh ..

    build="$prefix/build"
    mkdir -p "$build"

    cd "$build"

    rm -rfv "$prefix"/lib/python3

    if [ -n "$QLAT_MPICXX" ] ; then
        export MPICXX="$QLAT_MPICXX"
    fi
    export CXX="$MPICXX"
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

    prefix_python="$prefix/lib/python3"

    meson setup "$wd/qlat" \
        -Dpython.platlibdir="$prefix_python" -Dpython.purelibdir="$prefix_python" \
        --prefix="$prefix"

    time meson compile -j$num_proc
    meson install

    cd "$wd"

    mk-setenv.sh

    echo "!!!! $name build !!!!"

    rm -rf "$temp_dir" || true

} } 2>&1 | tee $prefix/log.$name.txt
