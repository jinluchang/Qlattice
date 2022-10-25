#!/bin/bash

. scripts/conf.sh

name=ninja

{

    time {

    echo "!!!! build $name !!!!"

    rm -rf $src_dir || true
    mkdir -p $src_dir || true
    cd $src_dir
    tar xaf $distfiles/$name-*

    cd $name-*

    python3 configure.py --bootstrap

    mkdir -p $prefix/bin
    cp -pv ninja $prefix/bin

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
