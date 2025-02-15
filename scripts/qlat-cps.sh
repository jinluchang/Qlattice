#!/usr/bin/env bash

name=qlat-cps

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    if [ -z "$(find-library.sh libcps.a)" ] ; then
        echo "qlat-cps: Cannot find 'libcps.a'. Stop."
        exit 1
    fi

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

    time-run meson setup "$wd/qlat-cps" \
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

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
    touch "$prefix"/build-successfully.txt
} } 2>&1 | tee $prefix/log.$name.txt
