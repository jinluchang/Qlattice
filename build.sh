#!/bin/bash

# This script would install LQPS in $lqps and all the dependencies in $prefix.
# Change $lqps if you need to install in other directories.
#
# CAUTION! This script could remove files in $lqps silently. Do not put anything
# important there. You have been warned.
#
# Authored by Luchang Jin

. conf.sh

if [ -z "$build_lqps" ] ; then
    build_lqps=true
fi
if [ -z "$build_libs" ] ; then
    build_libs=true
fi

if $build_libs ; then
    if [ -e $prefix ] ; then
        echo "$prefix already exist, continue to build will erase all its contents."
        echo "Use ./scripts/lqps.sh to build LQPS only."
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
    else
        ./scripts/gsl.sh
        ./scripts/fftw.sh
        ./scripts/cmake.sh
        ./scripts/eigen.sh
        ./scripts/hwloc.sh
        ./scripts/jansson.sh
        ./scripts/netloc.sh
        ./scripts/openmpi.sh
    fi
fi

if $build_lqps ; then
    ./scripts/timer.sh
    ./scripts/lqps.sh
    ./scripts/setenv.sh
fi

rm -rf $temp_dir
