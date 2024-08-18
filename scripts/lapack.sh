#!/usr/bin/env bash

name=lapack

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    rm -rf $src_dir
    mkdir -p $src_dir
    cd $src_dir
    time-run tar xaf $distfiles/$name-*

    rm -rf $build_dir || true
    mkdir -p $build_dir || true
    cd $build_dir

    time-run cmake $src_dir/$name-* \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DLAPACKE=ON \
        -DLAPACKE_WITH_TMG=ON \
        .

    time-run make -j$num_proc
    time-run make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
