#!/bin/bash

# This script would install Qlattice and all its dependencies in ``$prefix''.
# Change $prefix environment variable if you need to install in other directories.
#
# CAUTION! This script could remove files in $prefix silently. Do not put anything
# important there. You have been warned.
#
# Authored by Luchang Jin

echo "Need to run ./scripts/download.sh to download all needed packages."

. conf.sh

if which mpic++ >/dev/null 2>&1 ; then
    type mpic++
elif which mpicxx >/dev/null 2>&1 ; then
    type mpicxx
else
    echo "NO mpic++ or mpicxx available. Quit."
    exit
fi

if [ -e "$prefix" ] ; then
    echo "$prefix already exist, continue to build will erase all its contents."
    echo "Use ./scripts/qlat.sh to build Qlat only."
    echo "Ctrl-C to stop."
    for i in {10..0} ; do
        echo -n "$i "
        sleep 1;
    done
fi

rm -rf "$prefix" || true
mkdir -p "$prefix"

target="$1"

if [ -z "$target" ] ; then
    target=default
fi

if [ "$target" = default ] ; then
    ./scripts/setenv.sh
    . "$prefix"/setenv.sh
    #
    # ./scripts/gsl.sh
    # ./scripts/cmake.sh
    # ./scripts/lapack.sh
    ./scripts/fftw.sh
    ./scripts/fftwf.sh
    ./scripts/cuba.sh
    ./scripts/zlib.sh
    ./scripts/c-lime.sh
    ./scripts/eigen.sh
    ./scripts/qlat.sh
    #
    ./scripts/grid-avx2.sh
    ./scripts/gpt.sh
elif [ "$target" = sse4 ] ; then
    ./scripts/setenv.sh
    . "$prefix"/setenv.sh
    #
    ./scripts/fftw.sh
    ./scripts/fftwf.sh
    ./scripts/cuba.sh
    ./scripts/zlib.sh
    ./scripts/c-lime.sh
    ./scripts/eigen.sh
    ./scripts/qlat.sh
    #
    ./scripts/grid-sse4.sh
    ./scripts/gpt.sh
elif [ "$target" = gen16 ] ; then
    ./scripts/setenv.sh
    . "$prefix"/setenv.sh
    #
    ./scripts/fftw.sh
    ./scripts/fftwf.sh
    ./scripts/cuba.sh
    ./scripts/zlib.sh
    ./scripts/c-lime.sh
    ./scripts/eigen.sh
    ./scripts/qlat.sh
    #
    ./scripts/grid-gen-16.sh
    ./scripts/gpt.sh
elif [ "$target" = bnlknl ] ; then
    ./scripts/setenv-bnlknl.sh
    . "$prefix"/setenv.sh
    #
    ./scripts/fftw.sh
    ./scripts/fftwf.sh
    ./scripts/cuba.sh
    ./scripts/zlib.sh
    ./scripts/c-lime.sh
    ./scripts/eigen.sh
    ./scripts/qlat.sh
    #
    ./scripts/grid-knl.sh
    ./scripts/gpt.sh
elif [ "$target" = summit ] ; then
    ./scripts/setenv-summit.sh
    . "$prefix"/setenv.sh
    #
    ./scripts/fftw.sh
    ./scripts/fftwf.sh
    ./scripts/cuba.sh
    ./scripts/zlib.sh
    ./scripts/c-lime.sh
    ./scripts/eigen.sh
    ./scripts/qlat.sh
    #
    ./scripts/mpfr.sh
    ./scripts/grid-summit.sh
    ./scripts/gpt.sh
fi

rm -rf $temp_dir
