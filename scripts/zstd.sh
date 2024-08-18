#!/usr/bin/env bash

name=zstd

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    rm -rf $src_dir || true
    mkdir -p $src_dir || true
    cd $src_dir
    time-run tar xaf $distfiles/$name-*

    rm -rf $build_dir || true
    mkdir -p $build_dir || true

    cd $src_dir/$name-*/build/meson

    time-run meson setup $build_dir \
        -Dbin_programs=true -Dbin_contrib=true \
        --prefix=$prefix

    time-run meson compile -C $build_dir
    time-run meson install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
