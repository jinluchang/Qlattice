#!/bin/bash

name=zlib

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    rm -rf $src_dir || true
    mkdir -p $src_dir || true
    cd $src_dir
    debug tar xaf $distfiles/$name-*

    rm -rf $build_dir || true
    mkdir -p $build_dir || true
    cd $build_dir

    CFLAGS="-O3 -fPIC"

    debug $src_dir/$name-*/configure \
        --prefix=$prefix

    make -j$num_proc
    make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
