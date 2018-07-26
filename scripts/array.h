#!/bin/bash

. conf.sh

name=utils
echo "!!!! build $name !!!!"

mkdir -p $prefix/include/utils
cp -v $distfiles/compatible-headers/array.h $prefix/include/utils/array-compatible.h

echo "!!!! $name build !!!!"

rm -rf $temp_dir
