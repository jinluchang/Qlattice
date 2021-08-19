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

if [ -e $prefix ] ; then
    echo "$prefix already exist, continue to build will erase all its contents."
    echo "Use ./scripts/qlat.sh to build Qlat only."
    echo "Ctrl-C to stop."
    for i in {10..0} ; do
        echo -n "$i "
        sleep 1;
    done
fi

rm -rf $prefix/{bin,include,lib,pylib,share} || true
mkdir -p $prefix

if which mpic++ >/dev/null 2>&1 || which mpicxx >/dev/null 2>&1 ; then
    echo "mpi already exist. proceed."
else
    echo "NO mpic++ or mpicxx available. Quit."
    exit
fi

if [ -f $prefix/setenv.sh ] ; then
    echo "Setup file '$prefix/setenv.sh' already exist."
else
    echo "Install default '$prefix/setenv.sh' setup file."
    ./scripts/setenv.sh
fi


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

# ./scripts/grid-avx2.sh
# ./scripts/gpt.sh

rm -rf $temp_dir
