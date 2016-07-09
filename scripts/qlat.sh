#!/bin/bash

. conf.sh

name=qlat
echo "!!!! build $name !!!!"

mkdir -p $prefix/include
cp -v -r qlat $prefix/include

echo "!!!! $name build !!!!"

rm -rf $temp_dir
