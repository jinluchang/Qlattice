#!/bin/bash

name=bzip2

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    rm -rf $src_dir || true
    mkdir -p $src_dir || true
    cd $src_dir
    tar xaf $distfiles/$name-*

    cd $src_dir/$name-*
    make -j$num_proc CC="$CC -fPIC"
    make install PREFIX="$prefix"

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
