#!/bin/bash

name=c-lime

source qcore/set-prefix.sh $name

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh ..

    rm -rf $src_dir || true
    mkdir -p $src_dir || true
    cd $src_dir
    debug tar xzf $distfiles/$name.tar.gz

    cd *"$name"*
    ./autogen.sh
    cd ..

    rm -rf $build_dir || true
    mkdir -p $build_dir || true
    cd $build_dir

    export CFLAGS="$CFLAGS -fPIC"
    export CXXFLAGS="$CXXFLAGS -fPIC"

    debug "$src_dir"/*"$name"*/configure \
        --prefix="$prefix"

    make -j$num_proc
    make install

    cd "$wd"

    mk-setenv.sh

    echo "!!!! $name build !!!!"

    rm -rf "$temp_dir" || true

} } 2>&1 | tee $prefix/log.$name.txt
