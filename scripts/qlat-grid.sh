#!/usr/bin/env bash

name=qlat-grid

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    if ! which grid-config >/dev/null 2>&1 ; then
        echo "qlat-grid: Cannot find 'grid-config'. Stop."
        exit 1
    fi

    build="$prefix/build"
    mkdir -p "$build"

    cd "$build"

    CXX_ARR=($(grid-config --cxx))
    export CXX="${CXX_ARR[0]}"
    export CXXFLAGS="${CXX_ARR[@]:1} $CXXFLAGS"
    export LDFLAGS="${CXX_ARR[@]:1} $LDFLAGS"

    # export CXX_LD="$(grid-config --cxxld)"

    if [ -n "$QLAT_MPICXX" ] ; then
        export MPICXX="$QLAT_MPICXX"
        export CXX="$MPICXX"
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

    rm -rfv "$prefix"/build-successfully.txt

    echo "CXX=$CXX"
    echo "CXX_LD=$CXX_LD"
    echo "MPICXX=$MPICXX"
    echo "CXXFLAGS=$CXXFLAGS"
    echo "LDFLAGS=$LDFLAGS"
    echo "LIBS=$LIBS"

    is_fail=false

    time-run meson setup "$wd/qlat-grid" \
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
