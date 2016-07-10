#!/bin/bash

# This script would install QLAT in $qlat and all the dependencies in $prefix.
# Change $qlat if you need to install in other directories.
#
# CAUTION! This script could remove files in $qlat silently. Do not put anything
# important there. You have been warned.
#
# Authored by Luchang Jin

. conf.sh

if [ -z "$build_qlat" ] ; then
    build_qlat=true
fi
if [ -z "$build_libs" ] ; then
    build_libs=true
fi

if $build_libs ; then
    if [ -e $prefix ] ; then
        echo "$prefix already exist, continue to build will erase all its contents."
        echo "Use ./scripts/qlat.sh to build QLAT only."
        echo "Ctrl-C to stop."
        for i in {10..0} ; do
            echo -n "$i "
            sleep 1;
        done
    fi
	rm -rf $prefix || true
	mkdir -p $prefix
fi

if $build_libs ; then
    if [ $arch = bgq ] ; then
        ./scripts/gsl.sh
        ./scripts/fftw.sh
        ./scripts/cmake.sh
        ./scripts/lapack.sh
        ./scripts/eigen.sh
        ./scripts/hash-cpp.sh
    else
        ./scripts/gsl.sh
        ./scripts/fftw.sh
        ./scripts/cmake.sh
        ./scripts/eigen.sh
        ./scripts/hash-cpp.sh
        ./scripts/hwloc.sh
        ./scripts/jansson.sh
        ./scripts/netloc.sh
        ./scripts/openmpi.sh
    fi
fi

if $build_qlat ; then
    ./scripts/timer.sh
    ./scripts/qlat.sh
    ./scripts/setenv.sh
fi

rm -rf $temp_dir
