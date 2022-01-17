#!/bin/bash

. conf.sh

name=endian

{

echo "!!!! build $name !!!!"

mkdir -pv $prefix/include

cp qutils/compatible-endian.h $prefix/include/endian.h

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
