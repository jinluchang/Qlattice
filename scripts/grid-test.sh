#!/bin/bash

. scripts/conf.sh

name=grid-test

{

echo "!!!! build $name !!!!"

# grid_options="--dslash-asm --shm-hugepages --shm 4050"
# grid_options="--dslash-asm"
grid_options=""

geo_options="--grid 16.16.16.16 --mpi 1.1.1.1"

OMP_NUM_THREADS=4 $prefix/Grid/build/benchmarks/Benchmark_dwf_fp32 $grid_options $geo_options

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
