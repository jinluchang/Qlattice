#!/bin/bash

. conf.sh

name=omp

{

echo "!!!! build $name !!!!"

mkdir -pv $prefix/include

cp qutils/compatible-omp.h $prefix/include/omp.h

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
