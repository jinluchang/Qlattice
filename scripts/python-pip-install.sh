#!/bin/bash

name=python-pip-install

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    time-run pip3 install -U meson
    time-run pip3 install -U cython
    time-run pip3 install -U psutil
    time-run pip3 install -U sympy
    time-run pip3 install -U numpy
    time-run pip3 install -U scipy
    time-run pip3 install -U mpi4py

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
