#!/usr/bin/env bash

name=qlat

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    build="$prefix/build"
    mkdir -p "$build"

    cd "$build"

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

    rm -rfv "$prefix"/build-successfully.txt

    is_fail=false

    time-run meson setup "$wd/qlat" \
        --prefix="$prefix" \
        -Dpython.platlibdir="$prefix/lib/python3/qlat-packages" \
        -Dpython.purelibdir="$prefix/lib/python3/qlat-packages" \
        || is_fail=true

    if $is_fail ; then
        echo "cat meson-logs/meson-log.txt"
        cat meson-logs/meson-log.txt
        exit 1
    fi

    time-run meson compile -j$num_proc

    rm -rfv "$prefix"/bin
    rm -rfv "$prefix"/lib
    time-run meson install

    touch "$prefix"/build-successfully.txt

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
