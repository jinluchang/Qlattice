#!/bin/bash

. conf.sh

name=gpt-test

{

echo "!!!! build $name !!!!"

# grid_options="--dslash-asm --shm-hugepages --shm 4050"
# grid_options="--dslash-asm"
grid_options=""

geo_options="--grid 16.16.16.16 --mpi 1.1.1.1"

OMP_NUM_THREADS=4 $prefix/gpt/benchmarks/dslash.py $grid_options $geo_options --Ls 12 --N 10 

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
