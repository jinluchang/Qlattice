#!/bin/bash

. scripts/conf.sh

name=qlat

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/build-qlat"
    mkdir -p "$build"

    cd "$build"

    rm -rfv "$prefix"/include/qlat
    rm -rfv "$prefix"/include/qlat-setup.h

    meson "$wd/qlat" --prefix="$prefix"
    ninja -j$num_proc
    ninja install

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} |& tee $prefix/log.$name.txt
