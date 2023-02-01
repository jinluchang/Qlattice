#!/bin/bash

name=python-packages

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    if [ -z ${NPY_BLAS_ORDER+x} ] ; then
        export NPY_BLAS_ORDER=openblas
    fi

    if [ -z ${NPY_LAPACK_ORDER+x} ] ; then
        export NPY_LAPACK_ORDER=openblas
    fi

    if [ -z ${NPY_NUM_BUILD_JOBS+x} ] ; then
        export NPY_NUM_BUILD_JOBS=$num_proc
    fi

    find ~/.cache/pip/wheels -type f || true
    # rm -rfv ~/.cache/pip/wheels || true

    opts="--verbose --no-index --no-build-isolation --no-cache-dir -f $distfiles/python-packages"

    time-run pip3 install $opts psutil
    time-run pip3 install $opts sympy
    time-run pip3 install $opts cython
    time-run pip3 install $opts pythran
    time-run pip3 install $opts pybind11
    time-run pip3 install $opts numpy
    time-run pip3 install $opts scipy

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
