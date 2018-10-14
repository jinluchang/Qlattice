#!/bin/bash

. conf.sh

name=utils
echo "!!!! build $name !!!!"

mkdir -p $prefix/include/utils
cp -v utils/* $prefix/include/utils

echo "!!!! $name build !!!!"

rm -rf $temp_dir
