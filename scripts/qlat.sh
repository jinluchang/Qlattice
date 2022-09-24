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
    rm -rfv "$prefix"/lib/python3*/*-packages/cqlat.*
    rm -rfv "$prefix"/lib/python3*/*-packages/qlat
    rm -rfv "$prefix"/lib/python3*/*-packages/qlat_gpt.py
    rm -rfv "$prefix"/lib/python3*/*-packages/rbc_ukqcd*
    rm -rfv "$prefix"/lib/python3*/*-packages/auto_contractor*

    touch "$wd"/qlat/meson.build

    meson "$wd/qlat" --prefix="$prefix"
    ninja -j$num_proc
    ninja install

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} |& tee $prefix/log.$name.txt
