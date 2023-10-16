#!/bin/bash

name=python-jupyter

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    find ~/.cache/pip/wheels -type f || true
    # rm -rfv ~/.cache/pip/wheels || true

    opts="--verbose --upgrade"

    time-run pip3 install $opts pip || true
    time-run pip3 install $opts h5py || true
    time-run pip3 install $opts pandas || true
    time-run pip3 install $opts xarray || true
    time-run pip3 install $opts matplotlib || true
    time-run pip3 install $opts vpython || true
    time-run pip3 install $opts plotly || true
    time-run pip3 install $opts seaborn || true
    time-run pip3 install $opts notebook || true
    time-run pip3 install $opts dash || true
    time-run pip3 install $opts virtualenv || true
    time-run pip3 install $opts squarify || true
    time-run pip3 install $opts mathjax || true
    time-run pip3 install $opts lz4 || true
    time-run pip3 install $opts pycuda || true
    time-run pip3 install $opts pytools || true
    time-run pip3 install $opts Sphinx || true
    time-run pip3 install $opts myst-parser || true
    time-run pip3 install $opts torch torchvision torchaudio || true
    time-run pip3 install $opts ipywidgets || true
    time-run pip3 install $opts transformers || true
    time-run pip3 install $opts xformers || true
    time-run pip3 install $opts jupyterlab || true
    time-run pip3 install $opts jupyterlab-dash || true
    time-run pip3 install $opts jupyterhub || true
    # time-run pip3 install $opts jupyterlab-vpython

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
