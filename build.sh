#!/bin/bash

# This script would install Qlattice and all its dependencies in ``$prefix''.
# Change $prefix environment variable if you need to install in other directories.
#
# CAUTION! This script could remove files in $prefix silently. Do not put anything
# important there. You have been warned.
#
# Authored by Luchang Jin

. conf.sh

if [ -e $prefix ] ; then
    echo "$prefix already exist, continue to build will erase all its contents."
    echo "Use ./scripts/qlat.sh to build Qlat only."
    echo "Ctrl-C to stop."
    for i in {10..0} ; do
        echo -n "$i "
        sleep 1;
    done
fi
rm -rf $prefix || true
mkdir -p $prefix

if [ $arch = amd64 ] ; then
    ./scripts/hwloc.sh
    ./scripts/jansson.sh
    ./scripts/netloc.sh
    ./scripts/openmpi.sh
fi

./scripts/openssl.sh
./scripts/gsl.sh
./scripts/fftw.sh
./scripts/cmake.sh
./scripts/lapack.sh
./scripts/eigen.sh
./scripts/hash-cpp.sh
./scripts/utils.sh
./scripts/qlat.sh
./scripts/setenv.sh

rm -rf $temp_dir
