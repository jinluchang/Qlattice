#!/bin/bash

name=python-jupyter

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    find ~/.cache/pip/wheels -type f || true
    rm -rfv ~/.cache/pip/wheels || true

    opts="--verbose"

    time-run pip3 install $opts notebook
    time-run pip3 install $opts jupyterlab
    time-run pip3 install $opts jupyterhub
    time-run pip3 install $opts pandas
    time-run pip3 install $opts plotly
    time-run pip3 install $opts squarify
    time-run pip3 install $opts dash
    time-run pip3 install $opts seaborn
    time-run pip3 install $opts matplotlib
    time-run pip3 install $opts lz4

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
