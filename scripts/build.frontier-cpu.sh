#!/usr/bin/env bash

set -e

{
    ./scripts/setenv.frontier-cpu.sh
    ./scripts/qcore.sh

    ./scripts/gsl.sh
    ./scripts/cuba.sh
    # ./scripts/zlib.sh
    ./scripts/eigen.sh
    # ./scripts/perl.sh
    # ./scripts/openssl.sh
    ./scripts/libffi.sh
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
    ./scripts/gmp.sh
    ./scripts/mpfr.sh
    ./scripts/fftw.sh
    ./scripts/c-lime.sh

    ./scripts/qmp.sh
    ./scripts/qio.sh
    ./scripts/cps.sh

    ./scripts/grid-clehner.avx2.sh
    ./scripts/gpt.sh

    ./scripts/qlat-utils.sh
    ./scripts/qlat.sh
    ./scripts/qlat-cps.sh
    ./scripts/qlat-grid.sh

    date
} 2>&1 | tee $prefix/log.build.txt
