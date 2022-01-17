#!/bin/bash

. conf.sh

name=omp-remove

{

echo "!!!! build $name !!!!"

mkdir -pv $prefix/include

rm -v $prefix/include/endian.h

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
