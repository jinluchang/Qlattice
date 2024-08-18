#!/usr/bin/env bash

name=Cuba

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p $src_dir
    cd $src_dir
    time-run tar xzf $distfiles/$name-*.tar.*

    export CFLAGS="$CFLAGS -fPIC"
    export CXXFLAGS="$CXXFLAGS -fPIC"

    cd $name-*
    time-run ./configure \
        --prefix=$prefix

        # --build="$(arch)" \

    time-run make # do not support parallel build
    time-run make install
    time-run cp -rpv demo-* $prefix/share

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
