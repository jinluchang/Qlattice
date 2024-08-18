#!/usr/bin/env bash

set -e

{
    ./scripts/setenv.default.sh
    ./scripts/qcore.sh

    ./scripts/xz.sh
    ./scripts/zlib.sh
    ./scripts/bzip2.sh
    ./scripts/tar.sh
    ./scripts/gsl.sh
    ./scripts/gmp.sh
    ./scripts/mpfr.sh
    ./scripts/mpc.sh
    ./scripts/isl.sh
    ./scripts/gcc.sh
    ./scripts/bison.sh
    ./scripts/binutils.sh
    ./scripts/perl.sh
    ./scripts/openssl.sh
    ./scripts/cmake.sh
    ./scripts/libffi.sh
    ./scripts/tcl.sh
    ./scripts/tk.sh
    ./scripts/gnuplot.sh
    ./scripts/ncurses.sh
    ./scripts/libevent.sh
    ./scripts/tmux.sh

    ./scripts/openmpi.sh

    ./scripts/hdf5.sh
    ./scripts/openblas.sh
    ./scripts/fftw.sh
    ./scripts/cuba.sh
    ./scripts/eigen.sh

    rm -rfv ~/.cache/pip
    ./scripts/python.sh
    ./scripts/python-pip.sh
    ./scripts/re2c.sh
    ./scripts/ninja.sh
    ./scripts/ninja-script.sh
    ./scripts/python-meson.sh
    ./scripts/python-meson-py.sh
    ./scripts/python-packages.sh
    ./scripts/python-jupyter.sh

    ./scripts/autoconf.sh
    ./scripts/automake.sh
    ./scripts/c-lime.sh
    ./scripts/qmp.sh
    ./scripts/qio.sh
    ./scripts/cps.sh
    ./scripts/grid-clehner.avx2.sh
    ./scripts/gpt.sh

    ./scripts/qlat-packages.sh

    ./scripts/qlat-examples-py.sh
    ./scripts/qlat-examples-cpp.sh
    ./scripts/qlat-examples-py-gpt.sh
    ./scripts/qlat-examples-py-cps.sh
    ./scripts/qlat-examples-cpp-grid.sh

    ./scripts/hadrons.sh
    ./scripts/llvm-project.sh

    date
} 2>&1 | tee $prefix/log.build.txt
