#!/bin/bash

. scripts/res/conf.sh

name=endian

{

echo "!!!! build $name !!!!"

mkdir -pv $prefix/include

cp qlat-utils/include/qlat-utils/compatible-endian.h $prefix/include/endian.h

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name.txt
