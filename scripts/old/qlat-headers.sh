#!/usr/bin/env bash

name=qlat-headers

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    rm -rfv $prefix/include

    mkdir -pv $prefix/include
    cp -rpv qlat-utils/include/qlat-utils $prefix/include/
    cp -rpv qlat/include/qlat $prefix/include/
    cp -rpv qlat/include/qlat/qlat-setup.h $prefix/include/

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
