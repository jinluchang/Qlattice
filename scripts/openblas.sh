#!/usr/bin/env bash

name=OpenBLAS

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    rm -rf $src_dir || true
    mkdir -p $src_dir || true
    cd $src_dir
    time-run tar xaf $distfiles/$name-*

    cd $src_dir/$name-*

    time-run make MAKE_NB_JOBS=$num_proc DYNAMIC_ARCH=1
    time-run make install PREFIX=$prefix

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
