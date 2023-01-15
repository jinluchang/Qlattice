#!/bin/bash

. scripts/res/conf.sh

name=omp

{

echo "!!!! build $name !!!!"

mkdir -pv $prefix/include

cp qlat-utils/include/qlat-utils/compatible-omp.h /compatible-omp.h $prefix/include/omp.h

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name.txt
