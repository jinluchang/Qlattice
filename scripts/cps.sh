#!/usr/bin/env bash

name=CPS

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p "$prefix" || true
    time-run rsync -a --delete $distfiles/$name "$prefix"/

    export CXXFLAGS="$CXXFLAGS -fPIC -w -Wno-psabi"
    export CFLAGS="$CFLAGS -fPIC -w -Wno-psabi"

    opts=""
    if [ -n "$(find-library.sh libgmp.a)" ] ; then
        opts+=" --enable-gmp=$(find-library.sh libgmp.a)"
    fi
    if [ -n "$(find-library.sh libmpfr.a)" ] ; then
        opts+=" --enable-mpfr=$(find-library.sh libmpfr.a)"
    fi
    # if [ -n "$(find-library.sh libfftw3.a)" ] ; then
    #     opts+=" --enable-fftw=$(find-library.sh libfftw3.a)"
    # fi
    if [ -n "$(find-library.sh libqmp.a)" ] ; then
        opts+=" --enable-qmp=$(find-library.sh libqmp.a)"
    fi
    if [ -n "$(find-library.sh libqio.a)" ] ; then
        opts+=" --enable-qio=$(find-library.sh libqio.a)"
    fi

    if [ -n "$(find-library.sh libz.a)" ] ; then
        LDFLAGS+=" -L$(find-library.sh libz.a)/lib"
        CFLAGS+=" -I$(find-library.sh libz.a)/include"
        export LDFLAGS
        export CFLAGS
    fi

    export CC=$MPICC
    export CXX=$MPICXX

    prefix_actual="$(readlink -m "$prefix")"

    rm -rfv "$prefix_actual/build"
    mkdir -p "$prefix_actual/build"
    cd "$prefix_actual/build"

    time-run "$prefix_actual/$name"/cps_pp/configure \
        $opts \
        --prefix="$prefix_actual"

    time-run make -j$num_proc

    rm -rfv "$prefix/include"
    rm -rfv "$prefix/lib"
    mkdir -p "$prefix/include"
    mkdir -p "$prefix/lib"
    cp -rpv "$prefix/$name/cps_pp/include/"* "$prefix/include/"
    cp -rpv "$prefix/build/"*.h "$prefix/include/"
    cp -rpv "$prefix/build/include/"* "$prefix/include/"
    cp -rpv "$prefix/build/libcps.a" "$prefix/lib/"

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee "$prefix/log.$name.txt"
