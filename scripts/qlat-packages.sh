#!/bin/bash

. scripts/conf.sh

name=qlat-packages

{

    time {

    echo "!!!! build $name !!!!"

    find ~/.cache/pip/wheels -type f || true

    build="$prefix/qlat-packages"
    rm -rfv "$build" || true
    mkdir -p "$build"

    cd "$build"

    opts="--verbose --no-index --no-build-isolation --no-cache-dir -f $distfiles/python-packages"

    pip3 install $opts flit_core
    pip3 install $opts tomli
    pip3 install $opts wheel
    pip3 install $opts packaging
    pip3 install $opts pep517
    pip3 install $opts build
    pip3 install $opts pyproject-metadata
    pip3 install $opts setuptools_scm
    pip3 install $opts scikit-build
    pip3 install $opts meson-python
    pip3 install $opts psutil
    pip3 install $opts cython
    pip3 install $opts setuptools
    pip3 install $opts numpy
    pip3 install $opts sympy

    python3 -m build -ns -o "$build" "$wd"/qlat-utils

    pip3 install $opts -f "$build" qlat-utils

    python3 -m build -ns -o "$build" "$wd"/qlat

    pip3 install $opts -f "$build" qlat

    if grid-config --prefix >/dev/null 2>&1 ; then

        python3 -m build -ns -o "$build" "$wd"/qlat-grid

        pip3 install $opts -f "$build" qlat-grid

    fi

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
