#!/bin/bash

name=python-packages

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    if [ -z ${OPENBLAS+x} ] ; then
        export OPENBLAS="$(pkg-config openblas --variable=libdir)"
    fi

    if [ -z ${NPY_BLAS_ORDER+x} ] ; then
        export NPY_BLAS_ORDER=openblas
    fi

    if [ -z ${NPY_LAPACK_ORDER+x} ] ; then
        export NPY_LAPACK_ORDER=openblas
    fi

    if [ -z ${NPY_NUM_BUILD_JOBS+x} ] ; then
        export NPY_NUM_BUILD_JOBS=$num_proc
    fi

    if [ -z ${HDF5_DIR+x} ] ; then
        export HDF5_DIR="$(find-library.sh libhdf5)"
    fi

    find ~/.cache/pip/wheels -type f || true
    # rm -rfv ~/.cache/pip/wheels || true

    echo OPENBLAS="$OPENBLAS"
    echo NPY_BLAS_ORDER="$NPY_BLAS_ORDER"
    echo NPY_LAPACK_ORDER="$NPY_LAPACK_ORDER"
    echo NPY_NUM_BUILD_JOBS="$NPY_NUM_BUILD_JOBS"

    opts="--verbose --upgrade --no-index --no-build-isolation --no-cache-dir -f $distfiles/python-packages"

    time-run pip3 install $opts mpi4py
    time-run pip3 install $opts psutil
    time-run pip3 install $opts sympy
    time-run pip3 install $opts cython
    time-run pip3 install $opts pybind11
    time-run pip3 install $opts numpy
    time-run pip3 install $opts pythran
    time-run pip3 install $opts scipy
    time-run pip3 install $opts poetry_core
    time-run pip3 install $opts pkgconfig

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
