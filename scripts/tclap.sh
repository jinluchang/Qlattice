#!/bin/bash

name=tclap

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    rm -rf $src_dir || true
    mkdir -p $src_dir || true
    cd $src_dir
    debug tar xaf $distfiles/$name-*

    cd $src_dir/$name-*

    debug ./configure \
        --prefix=$prefix

    make -j$num_proc
    make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
