#!/usr/bin/env bash

name=python-jupyter

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    find ~/.cache/pip/wheels -type f || true
    # rm -rfv ~/.cache/pip/wheels || true

    opts="--verbose --upgrade"

    time-run pip3 install $opts pip || true

    opts="--verbose"

    time-run pip3 install $opts mpi4py || true
    time-run pip3 install $opts psutil || true
    time-run pip3 install $opts sympy || true
    time-run pip3 install $opts cython || true
    time-run pip3 install $opts pybind11 || true
    time-run pip3 install $opts numpy || true
    time-run pip3 install $opts pythran || true
    time-run pip3 install $opts scipy || true
    time-run pip3 install $opts poetry_core || true
    time-run pip3 install $opts pkgconfig || true

    time-run pip3 install $opts virtualenv || true
    time-run pip3 install $opts h5py || true
    time-run pip3 install $opts pandas || true
    time-run pip3 install $opts scikit-learn || true
    time-run pip3 install $opts xarray || true
    time-run pip3 install $opts matplotlib || true
    time-run pip3 install $opts plotly || true
    time-run pip3 install $opts seaborn || true
    time-run pip3 install $opts notebook || true
    time-run pip3 install $opts dash || true
    time-run pip3 install $opts squarify || true
    time-run pip3 install $opts mathjax || true
    time-run pip3 install $opts lz4 || true
    time-run pip3 install $opts pycuda || true
    time-run pip3 install $opts pytools || true
    time-run pip3 install $opts Sphinx || true
    time-run pip3 install $opts myst-parser || true
    time-run pip3 install $opts multiprocess dill || true
    time-run pip3 install $opts torch || true
    time-run pip3 install $opts jax jaxlib || true
    time-run pip3 install $opts torchvision torchaudio xformers || true
    time-run pip3 install $opts numba || true
    # time-run pip3 install $opts torch torchvision torchaudio xformers "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    time-run pip3 install $opts ipywidgets || true
    time-run pip3 install $opts transformers || true
    time-run pip3 install $opts jupyterlab || true
    time-run pip3 install $opts jupyterlab-dash || true
    time-run pip3 install $opts jupyterhub || true
    time-run pip3 install $opts vpython || true
    # time-run pip3 install $opts jupyterlab-vpython

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
