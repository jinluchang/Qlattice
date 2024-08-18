#!/usr/bin/env bash

set -e

{
    ./scripts/setenv.default.sh
    ./scripts/qcore.sh

    export CC=gcc
    export CXX=g++

    ./scripts/xz.sh
    ./scripts/tar.sh
    ./scripts/gmp.sh
    ./scripts/mpfr.sh
    ./scripts/mpc.sh
    ./scripts/gcc.sh
    ./scripts/binutils.sh
    ./scripts/cmake.sh
    ./scripts/perl.sh
    ./scripts/openblas.sh
    ./scripts/openssl.sh
    ./scripts/libffi.sh
    ./scripts/zlib.sh
    ./scripts/python.sh
    ./scripts/python-pip.sh
    ./scripts/python-packages.sh
    ./scripts/re2c.sh
    ./scripts/ninja.sh
    ./scripts/ninja-script.sh
    ./scripts/python-meson.sh
    ./scripts/python-meson-py.sh
    ./scripts/python-packages.sh

    ./scripts/llvm-project.sh

    export CC=
    export CXX=

    ./scripts/fftw.sh
    ./scripts/cuba.sh
    ./scripts/zlib.sh
    ./scripts/eigen.sh

    ./scripts/autoconf.sh
    ./scripts/automake.sh
    ./scripts/c-lime.sh
    ./scripts/hdf5.sh

    ./scripts/qlat-utils.sh
    ./scripts/qlat.sh

    # export CC=clang
    # export CXX=clang++
    # export CXX=mpiicpc

    ./scripts/grid-tblum.avx512.sh
    ./scripts/hadrons-tblum.sh

    date
} 2>&1 | tee $prefix/log.build.txt
