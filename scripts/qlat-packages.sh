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

    pip3 install $opts -f "$build" flit_core
    pip3 install $opts -f "$build" tomli
    pip3 install $opts -f "$build" wheel
    pip3 install $opts -f "$build" packaging
    pip3 install $opts -f "$build" pep517
    pip3 install $opts -f "$build" build
    pip3 install $opts -f "$build" pyproject-metadata
    pip3 install $opts -f "$build" setuptools_scm
    pip3 install $opts -f "$build" scikit-build
    pip3 install $opts -f "$build" meson-python
    pip3 install $opts -f "$build" psutil
    pip3 install $opts -f "$build" cython
    pip3 install $opts -f "$build" setuptools
    pip3 install $opts -f "$build" numpy
    pip3 install $opts -f "$build" sympy

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
