#!/bin/bash

. scripts/conf.sh

name=python-packages

{

    time {

    echo "!!!! build $name !!!!"

    find ~/.cache/pip/wheels -type f || true

    # rm -rfv ~/.cache/pip/wheels || true

    opts="--verbose --no-index --no-build-isolation --no-cache-dir -f $distfiles/python-packages"

    time pip3 install $opts psutil
    time pip3 install $opts sympy
    time pip3 install $opts cython
    time pip3 install $opts pythran
    time pip3 install $opts pybind11
    time pip3 install $opts numpy
    time pip3 install $opts scipy

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
