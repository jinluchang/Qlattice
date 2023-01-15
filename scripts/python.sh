#!/bin/bash

. scripts/res/conf.sh

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

    if [ -f "$prefix/lib64/libffi.a" ] ; then
        export LDFLAGS="-L$prefix/lib64"
        export LIBS="-lffi"
    elif [ -f "$prefix/lib/libffi.a" ] ; then
        export LDFLAGS="-L$prefix/lib"
        export LIBS="-lffi"
    fi

    $src_dir/$name-*/configure \
        --prefix=$prefix

    make -j$num_proc
    make install

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
