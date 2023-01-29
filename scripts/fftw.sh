#!/bin/bash

name=fftw

source qcore/set-prefix.sh $name

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh ..

    rm -rf $src_dir || true
    mkdir -p $src_dir || true
    cd $src_dir
    tar xzf $distfiles/$name-*.tar.gz

    cd $name-*
    export CFLAGS="$CFLAGS -fPIC"
    export CXXFLAGS="$CXXFLAGS -fPIC"

    ./configure \
        --prefix=$prefix \
        --enable-shared

    make -j$num_proc
    make install

    make clean

    ./configure \
        --prefix=$prefix \
        --enable-float \
        --enable-shared

    make -j$num_proc
    make install

    cd "$wd"

    mk-setenv.sh

    echo "!!!! $name build !!!!"

    rm -rf "$temp_dir" || true

} } 2>&1 | tee "$prefix/log.$name.txt"
