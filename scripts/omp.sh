#!/bin/bash

. conf.sh

name=utils
echo "!!!! build $name !!!!"

mkdir -p $prefix/include/utils
cp -v $distfiles/compatible-headers/omp.h $prefix/include/utils

echo "!!!! $name build !!!!"

rm -rf $temp_dir
