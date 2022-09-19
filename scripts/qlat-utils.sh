#!/bin/bash

. scripts/conf.sh

name=qlat-utils

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/build-qlat-utils"
    mkdir -p "$build"

    cd "$build"

    meson "$wd/qlat-utils" --prefix="$prefix"
    time ninja -j$num_proc
    ninja install

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} |& tee $prefix/log.$name.txt
