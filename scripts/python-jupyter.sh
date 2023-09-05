#!/bin/bash

name=python-jupyter

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    find ~/.cache/pip/wheels -type f || true
    # rm -rfv ~/.cache/pip/wheels || true

    opts="--verbose"

    time-run pip3 install $opts h5py
    time-run pip3 install $opts pandas
    time-run pip3 install $opts xarray
    time-run pip3 install $opts matplotlib
    time-run pip3 install $opts vpython
    time-run pip3 install $opts plotly
    time-run pip3 install $opts seaborn
    time-run pip3 install $opts notebook
    time-run pip3 install $opts dash
    time-run pip3 install $opts virtualenv
    time-run pip3 install $opts squarify
    time-run pip3 install $opts mathjax
    time-run pip3 install $opts lz4
    time-run pip3 install $opts pycuda
    time-run pip3 install $opts pytools
    time-run pip3 install $opts Sphinx
    time-run pip3 install $opts myst-parser
    time-run pip3 install $opts torch torchvision torchaudio
    time-run pip3 install $opts ipywidgets
    time-run pip3 install $opts transformers
    time-run pip3 install $opts xformers
    time-run pip3 install $opts jupyterlab
    time-run pip3 install $opts jupyterlab-vpython
    time-run pip3 install $opts jupyterlab-dash
    time-run pip3 install $opts jupyterhub

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
