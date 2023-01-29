#!/bin/bash

name=grid-test

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    cd "$prefix"

    # grid_options="--dslash-asm --shm-hugepages --shm 4050"
    # grid_options="--dslash-asm"
    grid_options=""

    geo_options="--grid 16.16.16.16 --mpi 1.1.1.1"

    OMP_NUM_THREADS=4 $prefix/../Grid-clehner/bin/Benchmark_dwf_fp32 $grid_options $geo_options

    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
