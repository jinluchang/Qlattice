#!/usr/bin/env bash

name=python-pip

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    find ~/.cache/pip/wheels -type f || true
    # rm -rfv ~/.cache/pip/wheels || true

    opts="--verbose --upgrade --no-index --no-build-isolation --no-cache-dir -f $distfiles/python-packages"

    time-run pip3 install $opts flit_core
    time-run pip3 install $opts wheel
    time-run pip3 install $opts setuptools
    time-run pip3 install $opts pip

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
