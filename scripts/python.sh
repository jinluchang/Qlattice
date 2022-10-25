#!/bin/bash

. scripts/conf.sh

name=Python

{

    time {

    echo "!!!! build $name !!!!"

    rm -rf $src_dir || true
    mkdir -p $src_dir || true
    cd $src_dir
    tar xaf $distfiles/$name-*

    rm -rf $build_dir || true
    mkdir -p $build_dir || true
    cd $build_dir

    export LDFLAGS="-L$prefix/lib64 -L$prefix/lib"
    export LIBS="-lffi"

    $src_dir/$name-*/configure \
        --prefix=$prefix \
        --with-openssl=$prefix

    make -j$num_proc
    make install

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
