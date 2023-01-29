#!/bin/bash

name=python-jupyter

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    find ~/.cache/pip/wheels -type f || true

    rm -rfv ~/.cache/pip/wheels || true

    time pip3 install notebook
    time pip3 install jupyterlab
    time pip3 install jupyterhub
    time pip3 install pandas
    time pip3 install plotly
    time pip3 install squarify
    time pip3 install dash
    time pip3 install seaborn
    time pip3 install matplotlib
    time pip3 install lz4

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
