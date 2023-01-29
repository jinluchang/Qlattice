#!/bin/bash

name=python-jupyter

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    find ~/.cache/pip/wheels -type f || true
    rm -rfv ~/.cache/pip/wheels || true

    opts="--verbose"

    time pip3 install $opts notebook
    time pip3 install $opts jupyterlab
    time pip3 install $opts jupyterhub
    time pip3 install $opts pandas
    time pip3 install $opts plotly
    time pip3 install $opts squarify
    time pip3 install $opts dash
    time pip3 install $opts seaborn
    time pip3 install $opts matplotlib
    time pip3 install $opts lz4

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
