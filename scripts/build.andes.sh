#!/usr/bin/env bash

set -e

{
    ./scripts/setenv.andes.sh
    ./scripts/qcore.sh

    # ./scripts/gmp.sh
    # ./scripts/mpfr.sh
    # ./scripts/mpc.sh
    # ./scripts/isl.sh
    # ./scripts/gcc.sh
    # ./scripts/bison.sh
    # ./scripts/binutils.sh

    ./scripts/gsl.sh
    ./scripts/cuba.sh
    # ./scripts/zlib.sh
    ./scripts/eigen.sh
    # ./scripts/perl.sh
    ./scripts/openssl.sh
    ./scripts/libffi.sh
    ./scripts/ncurses.sh
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

    ./scripts/autoconf.sh
    ./scripts/automake.sh
    ./scripts/fftw.sh
    ./scripts/c-lime.sh

    ./scripts/qmp.sh
    ./scripts/qio.sh
    ./scripts/cps.sh

    ./scripts/grid-clehner.avx2.sh
    ./scripts/gpt.sh

    ./scripts/qlat-utils.sh
    ./scripts/qlat.sh
    ./scripts/qlat-grid.sh
    ./scripts/qlat-cps.sh

    date
} 2>&1 | tee $prefix/log.build.txt
