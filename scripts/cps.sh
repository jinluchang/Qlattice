#!/usr/bin/env bash

# Fail the build step (and CI) when the build fails: the trailing
# `2>&1 | tee` pipeline below would otherwise mask the exit status, which is
# how a gcc-15 CPS build failure went unnoticed while CI happily ran the CPS
# tests without a qlat_cps installation.
set -e -o pipefail

name=CPS

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p "$prefix" || true
    time-run rsync -a --delete $distfiles/$name "$prefix"/

    export CXXFLAGS="$CXXFLAGS -fPIC -w -Wno-psabi"
    # CPS is legacy C: gcc 15 turned some C diagnostics that used to be
    # warnings into errors by default in all C modes (e.g.
    # cps_pp/src/util/dirac_op/d_op_base/comsrc/matrix_dgemm.c:9,
    # `const double **pA = A;` hits -Wincompatible-pointer-types).  Demote just
    # those diagnostics back to warnings instead of weakening qlat's own
    # compiler flags; this is the same set nixpkgs/cps.nix already passes via
    # NIX_CFLAGS_COMPILE to build this exact CPS revision with gcc 15.  The -w
    # above keeps the CPS build log quiet.
    export CFLAGS="$CFLAGS -fPIC -w -Wno-psabi"
    export CFLAGS="$CFLAGS -Wno-error=incompatible-pointer-types"
    export CFLAGS="$CFLAGS -Wno-error=int-conversion"
    export CFLAGS="$CFLAGS -Wno-error=implicit-function-declaration"

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
