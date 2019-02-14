#!/bin/bash

. conf.sh

name=qutils
echo "!!!! build $name !!!!"

rm -rf $prefix/include/$name
mkdir -p $prefix/include/$name
cp -pv qutils/*.h $prefix/include/$name

echo "!!!! $name built !!!!"

rm -rf $temp_dir
