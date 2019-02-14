#!/bin/bash

. conf.sh

name=qlat
echo "!!!! build $name !!!!"

rm -rf $prefix/include/$name
mkdir -p $prefix/include/$name
cp -pv qlat/*.h $prefix/include/$name

echo "!!!! $name build !!!!"

rm -rf $temp_dir
