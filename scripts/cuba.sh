#!/bin/bash

name=Cuba

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p $src_dir
    cd $src_dir
    tar xzf $distfiles/$name-*.tar.*

    export CFLAGS="$CFLAGS -fPIC"
    export CXXFLAGS="$CXXFLAGS -fPIC"

    cd $name-*
    debug ./configure \
        --build="$(arch)" \
        --prefix=$prefix

    make # do not support parallel build
    make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
