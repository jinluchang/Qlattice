#!/bin/bash

. scripts/res/conf.sh

name=gpt-test

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/test-gpt"
    mkdir -p "$build"
    cd "$build"

    # grid_options="--dslash-asm --shm-hugepages --shm 4050"
    # grid_options="--dslash-asm"
    grid_options=""

    geo_options="--grid 16.16.16.16 --mpi 1.1.1.1"

    echo 'import numpy' | cat - $prefix/gpt/benchmarks/dslash.py > dslash.py
    OMP_NUM_THREADS=4 python3 dslash.py $grid_options $geo_options --Ls 12 --N 10

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
