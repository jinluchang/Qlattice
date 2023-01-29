#!/bin/bash

name=hdf5

source qcore/set-prefix.sh $name

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh ..

    rm -rf $src_dir || true
    mkdir -p $src_dir || true
    cd $src_dir
    tar xjf $distfiles/$name-*

    rm -rf $build_dir || true
    mkdir -p $build_dir || true
    cd $build_dir

    $src_dir/$name-*/configure \
        --prefix=$prefix \
        --enable-cxx

    make -j$num_proc
    make install

    cd "$wd"

    mk-setenv.sh

    echo "!!!! $name build !!!!"

    rm -rf "$temp_dir" || true

} } 2>&1 | tee $prefix/log.$name.txt
