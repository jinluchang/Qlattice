#!/bin/bash

. conf.sh

name=lqps
echo "!!!! build $name !!!!"

mkdir -p $prefix/include
cp -v -r lqps $prefix/include

echo "!!!! $name build !!!!"

rm -rf $temp_dir
