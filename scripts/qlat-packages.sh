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

    pip3 install $opts -f "$build" meson-python

    python3 -m build -ns -o "$build" "$wd"/qlat-utils

    pip3 install $opts -f "$build" qlat-utils

    python3 -m build -ns -o "$build" "$wd"/qlat

    pip3 install $opts -f "$build" qlat

    python3 -m build -ns -o "$build" "$wd"/qlat-grid

    pip3 install $opts -f "$build" qlat-grid

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
