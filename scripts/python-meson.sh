#!/usr/bin/env bash

name=python-meson

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    find ~/.cache/pip/wheels -type f || true
    # rm -rfv ~/.cache/pip/wheels || true

    opts="--verbose --upgrade --no-index --no-build-isolation --no-cache-dir -f $distfiles/python-packages"

    time-run pip3 install $opts meson
    time-run pip3 install $opts ninja

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
