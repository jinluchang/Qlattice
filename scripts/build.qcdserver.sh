#!/bin/bash

set -e

{

./scripts/setenv.default.sh
./scripts/qcore.sh

./scripts/xz.sh
./scripts/tar.sh
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
./scripts/zlib.sh
./scripts/bzip2.sh
./scripts/tcl.sh
./scripts/tk.sh
./scripts/openblas.sh
./scripts/python.sh
./scripts/python-pip.sh
./scripts/re2c.sh
./scripts/ninja.sh
./scripts/ninja-script.sh
./scripts/python-meson.sh
./scripts/python-meson-py.sh
./scripts/python-packages.sh

./scripts/openmpi.sh

./scripts/fftw_mpi.sh
./scripts/cuba.sh
./scripts/eigen.sh
./scripts/qlat.sh

./scripts/qlat-examples-py.sh
./scripts/qlat-examples-cpp.sh

./scripts/autoconf.sh
./scripts/automake.sh
./scripts/c-lime.sh
./scripts/hdf5.sh
./scripts/grid-clehner.avx2.sh
./scripts/gpt.sh
./scripts/qlat-grid.sh

./scripts/qlat-examples-py-gpt.sh

./scripts/gnuplot.sh
./scripts/python-jupyter.sh

./scripts/cmake.sh
./scripts/llvm-project.sh

} 2>&1 | tee $prefix/log.build.txt
