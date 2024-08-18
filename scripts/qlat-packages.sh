#!/usr/bin/env bash

name=qlat-packages

./scripts/qlat-packages-uninstall.sh

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    time-run find ~/.cache/pip/wheels -type f || true
    time-run find ~/.cache/pip/wheels -type f -name 'qlat*.whl' -delete || true

    build="$prefix/packages"
    time-run rm -rfv "$build" || true
    mkdir -p "$build"

    cd "$build"

    # opts="--verbose --no-index --no-cache-dir -f $distfiles/python-packages -f $build"
    # opts="--user --verbose --force-reinstall -f $build"
    opts="--verbose -f $build --force-reinstall --config-settings compile-args=-j$num_proc"

    time-run python3 -m build -ns -o "$build" "$wd"/qlat-utils
    time-run pip3 install $opts qlat-utils
    time-run python3 -m build -ns -o "$build" "$wd"/qlat
    time-run pip3 install $opts qlat
    if [ -n "$(find-library.sh libcps.a)" ] ; then
        time-run python3 -m build -ns -o "$build" "$wd"/qlat-cps
        time-run pip3 install $opts qlat-cps
    fi
    if grid-config --prefix >/dev/null 2>&1 ; then
        time-run python3 -m build -ns -o "$build" "$wd"/qlat-grid
        time-run pip3 install $opts qlat-grid
    fi

    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
