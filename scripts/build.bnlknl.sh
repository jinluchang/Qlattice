#!/usr/bin/env bash

set -e

{
    ./scripts/setenv.bnlknl.sh
    ./scripts/qcore.sh

    ./scripts/tar.sh
    ./scripts/zlib.sh
    ./scripts/xz.sh
    ./scripts/gsl.sh
    ./scripts/gmp.sh
    ./scripts/mpfr.sh
    ./scripts/mpc.sh
    ./scripts/isl.sh
    ./scripts/gcc.sh
    ./scripts/binutils.sh
    ./scripts/perl.sh
    ./scripts/openssl.sh
    ./scripts/libffi.sh
    ./scripts/openblas.sh
    ./scripts/hdf5.sh
    ./scripts/python.sh
    ./scripts/python-pip.sh
    ./scripts/re2c.sh
    ./scripts/ninja.sh
    ./scripts/ninja-script.sh
    ./scripts/python-meson.sh
    ./scripts/python-meson-py.sh
    ./scripts/python-packages.sh

    ./scripts/cuba.sh
    ./scripts/eigen.sh

    ./scripts/fftw.sh

    ./scripts/qlat-utils.sh
    ./scripts/qlat.sh

    ./scripts/autoconf.sh
    ./scripts/automake.sh

    ./scripts/c-lime.sh
    ./scripts/grid-clehner.knl.sh
    ./scripts/gpt.sh
    ./scripts/qlat-grid.sh

    ./scripts/qmp.sh
    ./scripts/qio.sh
    ./scripts/cps.sh
    ./scripts/qlat-cps.sh

    ./scripts/qlat-examples-py.sh
    ./scripts/qlat-examples-cpp.sh
    ./scripts/qlat-examples-py-gpt.sh
    ./scripts/qlat-examples-py-cps.sh
    ./scripts/qlat-examples-cpp-grid.sh

    date
} 2>&1 | tee $prefix/log.build.txt
