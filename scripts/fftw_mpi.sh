#!/usr/bin/env bash

name=fftw_mpi
pkgname=fftw

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    rm -rf $src_dir || true
    mkdir -p $src_dir || true
    cd $src_dir
    time-run tar xzf $distfiles/$pkgname-*.tar.gz

    cd $pkgname-*

    export CFLAGS="$CFLAGS -fPIC"
    export CXXFLAGS="$CXXFLAGS -fPIC"

    time-run ./configure \
        --prefix=$prefix \
        --enable-mpi \
        --enable-shared

    time-run make -j$num_proc
    time-run make install

    time-run make clean

    time-run ./configure \
        --prefix=$prefix \
        --enable-float \
        --enable-mpi \
        --enable-shared

    time-run make -j$num_proc
    time-run make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
