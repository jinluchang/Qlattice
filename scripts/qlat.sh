#!/bin/bash

. scripts/conf.sh

name=qlat

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/build-qlat"
    mkdir -p "$build"

    cd "$build"

    meson "$wd/qlat" --prefix="$prefix"
    ninja -j$num_proc
    ninja install

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} |& tee $prefix/log.$name.txt
