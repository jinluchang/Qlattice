#!/bin/bash

name=qlat-packages

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    time-run find ~/.cache/pip/wheels -type f || true
    # time-run rm -rfv ~/.cache/pip/wheels || true

    build="$prefix/packages"
    time-run rm -rfv "$build" || true
    mkdir -p "$build"

    cd "$build"

    # opts="--verbose --no-index --no-cache-dir -f $distfiles/python-packages -f $build"
    # opts="--user --verbose --force-reinstall -f $build"
    opts="--verbose -f $build"

    time-run pip3 uninstall -y qlat-utils qlat qlat-grid qlat-cps || true

    time-run python3 -m build -ns -o "$build" "$wd"/qlat-utils
    time-run pip3 install $opts qlat-utils
    time-run python3 -m build -ns -o "$build" "$wd"/qlat
    time-run pip3 install $opts qlat
    if grid-config --prefix >/dev/null 2>&1 ; then
        time-run python3 -m build -ns -o "$build" "$wd"/qlat-grid
        time-run pip3 install $opts qlat-grid
    fi
    if [ -n "$(find-library.sh libcps.a)" ] ; then
        time-run python3 -m build -ns -o "$build" "$wd"/qlat-cps
        time-run pip3 install $opts qlat-cps
    fi

    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
