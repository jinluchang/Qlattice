#!/usr/bin/env bash

name=qlat-test

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    cd "$prefix"

    OMP_NUM_THREADS=4 time-run $wd/examples-py/hmc-pure-gauge.py

    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
