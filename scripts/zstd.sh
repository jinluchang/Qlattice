#!/bin/bash

name=zstd

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

    cd $src_dir/$name-*/build/meson

    meson setup $build_dir \
        -Dbin_programs=true -Dbin_contrib=true \
        --prefix=$prefix

    meson compile -C $build_dir
    meson install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
