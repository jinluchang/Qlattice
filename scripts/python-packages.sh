#!/usr/bin/env bash

name=python-packages

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    find ~/.cache/pip/wheels -type f || true
    # rm -rfv ~/.cache/pip/wheels || true

    opts="--verbose --upgrade --no-index --no-build-isolation --no-cache-dir -f $distfiles/python-packages"

    time-run pip3 install $opts psutil || true
    time-run pip3 install $opts sympy || true
    time-run pip3 install $opts cython || true
    time-run pip3 install $opts mpi4py || true
    time-run pip3 install $opts pybind11 || true
    time-run pip3 install $opts numpy || true
    time-run pip3 install $opts pythran || true
    time-run pip3 install $opts scipy || true
    time-run pip3 install $opts poetry_core || true
    time-run pip3 install $opts pkgconfig || true

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
