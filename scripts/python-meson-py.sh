#!/bin/bash

name=python-meson-py

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    find ~/.cache/pip/wheels -type f || true
    rm -rfv ~/.cache/pip/wheels || true

    opts="--verbose --no-index --no-build-isolation --no-cache-dir -f $distfiles/python-packages"

    time pip3 install $opts flit_core
    time pip3 install $opts tomli
    time pip3 install $opts packaging
    time pip3 install $opts pep517
    time pip3 install $opts build
    time pip3 install $opts pyproject-metadata
    time pip3 install $opts setuptools_scm
    time pip3 install $opts scikit-build
    time pip3 install $opts meson-python

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
