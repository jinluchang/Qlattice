#!/bin/bash

. scripts/conf.sh

name=qlat-test

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/test-qlat"
    mkdir -p "$build"
    cd "$build"

    OMP_NUM_THREADS=4 $wd/examples-py/hmc-pure-gauge.py

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} |& tee $prefix/log.$name.txt
