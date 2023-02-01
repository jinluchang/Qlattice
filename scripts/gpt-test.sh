#!/bin/bash

name=gpt-test

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    cd "$prefix"

    # grid_options="--dslash-asm --shm-hugepages --shm 4050"
    # grid_options="--dslash-asm"
    grid_options=""

    geo_options="--grid 16.16.16.16 --mpi 1.1.1.1"

    echo 'import numpy' | cat - "$prefix"/../gpt/src/benchmarks/dslash.py > dslash.py
    OMP_NUM_THREADS=4 time-run python3 dslash.py $grid_options $geo_options --Ls 12 --N 10

    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
